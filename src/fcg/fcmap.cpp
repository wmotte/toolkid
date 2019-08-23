#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "rsCommon.h"
#include "tkdCmdParser.h"
#include <map>
#include "itkImageRegionConstIteratorWithIndex.h"
#include "rsPCA.h"
#include "itkImageRegionIterator.h"

int main( int argc, char ** argv )
{
  tkd::CmdParser p( "fcmap", "Calculate functional connectivity from seed ROI's" );

  std::string inputImageName;
  std::string labelImageName;
  std::string outputFileName;
  std::vector< int > labelRange;
  std::vector< int > ignore;
  bool rToZ = false;
  bool sampleCoG = false;
  bool doPCA = false;

  p.AddArgument( inputImageName, "input" )
    ->AddAlias( "i" )
    ->SetDescription( "Input 4D time series image" )
    ->SetRequired( true );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
    ->SetDescription( "Output connectivity maps (4D)" )
    ->SetRequired( true );

  p.AddArgument( labelImageName, "labels" )
    ->AddAlias( "l" )
    ->SetDescription( "Input seed ROI's in 3D label image" )
    ->SetRequired( true );

  p.AddArgument( labelRange, "label-range" )
    ->AddAlias( "lr" )
    ->SetDescription( "Specify a label range" )
    ->SetMinMax( 2, 2 );

  p.AddArgument( ignore, "ignore" )
    ->AddAlias( "ig" )
    ->SetDescription( "Ignore samples at begin and end" )
    ->SetMinMax( 2, 2 );

  p.AddArgument( rToZ, "fisher-z" )
    ->AddAlias( "z" )
    ->SetDescription( "Output Fisher z'-transformed connectivity maps" );

  p.AddArgument( sampleCoG, "cog" )
    ->AddAlias( "c" )
    ->SetDescription( "Seed vector from center-of-gravity" );

  p.AddArgument( doPCA, "pca" )
    ->SetDescription( "Extract characteristic ROI profiles with PCA" );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  typedef float PixelType;
  typedef itk::Image< PixelType, 4 > InputType;
  typedef itk::Image< PixelType, 4 > OutputType;
  typedef itk::Image< int, 3 > LabelType;

  typedef itk::ImageFileReader< InputType > InputReaderType;
  typedef itk::ImageFileReader< LabelType > LabelReaderType;
  typedef itk::ImageFileWriter< OutputType > WriterType;

  InputReaderType::Pointer inputReader = InputReaderType::New();
  inputReader->SetFileName( inputImageName.c_str() );
  inputReader->Update();

  InputType::Pointer input = inputReader->GetOutput();
  inputReader = 0;

  LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName( labelImageName.c_str() );
  labelReader->Update();

  LabelType::Pointer label = labelReader->GetOutput();
  labelReader = 0;

  std::map< int, int > labels;

  if ( labelRange.size() == 2 && labelRange[ 1 ] >= labelRange[ 0 ] )
    {
    int count = 0;
    for( int i = labelRange[ 0 ]; i <= labelRange[ 1 ]; ++i )
      {
      labels[ i ] = count++;
      }
    }
  else
    {
    int count = 0;
    itk::ImageRegionIterator< LabelType > it( label, label->GetLargestPossibleRegion() );
    while( !it.IsAtEnd() )
      {
      if ( it.Value() > 0 && labels.find( it.Value() ) == labels.end() )
        {
        labels[ it.Value() ] = count++;
        }
      ++it;
      }
    }

  int numberOfLabels = labels.size();

  InputType::RegionType region = input->GetLargestPossibleRegion();
  InputType::SizeType size = region.GetSize();
  int voxels = size[ 0 ] * size[ 1 ] * size[ 2 ];
  int timepoints = size[ 3 ];

  OutputType::SizeType outputSize = size;
  outputSize[ 3 ] = numberOfLabels;
  OutputType::RegionType outputRegion;
  outputRegion.SetSize( outputSize );


  OutputType::Pointer output = OutputType::New();
  output->CopyInformation( input );
  output->SetRegions( outputRegion );
  output->Allocate();
  output->FillBuffer( 0 );

  vnl_matrix_ref< PixelType > matrix( timepoints, voxels, input->GetPixelContainer()->GetBufferPointer() );
  const int* mask = label->GetPixelContainer()->GetBufferPointer();
  vnl_matrix_ref< PixelType > fc( numberOfLabels, voxels, output->GetPixelContainer()->GetBufferPointer() );

  vnl_matrix< PixelType > mean( timepoints, numberOfLabels );
  mean.fill( 0 );

  if ( sampleCoG )
    {
    std::vector< LabelType::PointType > centers;
    for( int i = 0; i < numberOfLabels; ++i )
      {
      LabelType::PointType point;
      for( int d = 0; d < 3; ++d )
        {
        point[ d ] = 0;
        }
      centers.push_back( point );
      }

    vnl_vector< int > count( numberOfLabels );
    count.fill( 0 );

    itk::ImageRegionConstIteratorWithIndex< LabelType > it( label, label->GetLargestPossibleRegion() );
    const int* firstPixel = label->GetPixelContainer()->GetBufferPointer();

    for( ; !it.IsAtEnd(); ++it )
      {
      int labelId = it.Get();

      if ( labelId < 1 )
        {
        continue;
        }

      labelId = labels[ labelId ];

      LabelType::IndexType index = it.GetIndex();
      LabelType::PointType point;
      label->TransformIndexToPhysicalPoint( index, point );

      for( int d = 0; d < 3; ++d )
        {
        centers[ labelId ][ d ] += point[ d ];
        }

      count[ labelId ] += 1;
      }

    for( int i = 0; i < numberOfLabels; ++i )
      {
      if ( count[ i ] == 0 )
        {
        continue;
        }

      for( int d = 0; d < 3; ++d )
        {
        centers[ i ][ d ] /= count[ i ];
        }

      LabelType::IndexType index;
      label->TransformPhysicalPointToIndex( centers[ i ], index );

      if ( !label->GetLargestPossibleRegion().IsInside( index ) )
        {
        continue;
        }

      std::cout << "Label " << i << " has center " << centers[ i ] << " and index " << index << std::endl;
      InputType::IndexType inputIndex;
      inputIndex[ 0 ] = index[ 0 ];
      inputIndex[ 1 ] = index[ 1 ];
      inputIndex[ 2 ] = index[ 2 ];
      for( int j = 0; j < timepoints; ++j )
        {
        inputIndex[ 3 ] = j;
        mean( j, i ) = input->GetPixel( inputIndex );
        }
      }
    }
  else if ( doPCA )
    {
    for( std::map< int, int >::iterator i = labels.begin(); i != labels.end(); ++i )
      {
      int label = i->first;
      int column = i->second;

      int numberOfSamples = 0;
      for( int voxel = 0; voxel < voxels; ++voxel )
        {
        if ( mask[ voxel ] != label )
          {
          continue;
          }

        ++numberOfSamples;
        }

      vnl_matrix< PixelType > sampleMatrix( timepoints, numberOfSamples );

      int sample = 0;
      for( int voxel = 0; voxel < voxels; ++voxel )
        {
        if ( mask[ voxel ] != label )
          {
          continue;
          }

        sampleMatrix.set_column( sample++, matrix.get_column( voxel ) );
        }

      std::cout << "Label " << label << ", " << sampleMatrix.rows() << "x" << sampleMatrix.cols() << std::endl;
      rs::PCA pca( sampleMatrix );
      pca.TransformPCA();
      vnl_vector< PixelType > signal = pca.GetComponent( 0 );
      mean.set_column( column, signal );
      }
    }
  else
    {
    vnl_vector< int > count( numberOfLabels );
    count.fill( 0 );
    for( int i = 0; i < voxels; ++i )
      {
      if ( mask[ i ] < 1 )
        {
        continue;
        }

      int column = labels[ mask[ i ] ];

      mean.set_column( column, mean.get_column( column ) + matrix.get_column( i ) );
      ++count[ column ];
      }

    for( int i = 0; i < numberOfLabels; ++i )
      {
      mean.set_column( i, mean.get_column( i ) / static_cast< PixelType >( count( i ) ) );
      }
    }

  int vectorLength = timepoints;
  int vectorStart = 0;

  if ( ignore.size() == 2 )
    {
    vectorStart = ignore[ 0 ];
    vectorLength = timepoints - vectorStart - ignore[ 1 ];
    }

  mean = mean.extract( vectorLength, timepoints, vectorStart, 0 );

  for( int i = 0; i < voxels; ++i )
    {
    const vnl_vector< PixelType > v = matrix.get_column( i ).extract( vectorLength, vectorStart );

    for( int j = 0; j < numberOfLabels; ++j )
      {
      PixelType cc = rs::Common< PixelType >::CorrelationCoefficient( v, mean.get_column( j ) );

      if ( cc != cc )
        {
        continue;
        }

      if ( rToZ )
        {
        cc = rs::Common< PixelType >::RtoZ( cc );
        }

      fc( j, i ) = cc;
      }
    }

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputFileName.c_str() );
  writer->SetInput( output );
  writer->Update();

  return 0;
}
