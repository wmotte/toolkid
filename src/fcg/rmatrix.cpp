#include "rmatrix.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include <sstream>
#include "tkdCmdParser.h"

namespace fcg
{

void RMatrix::Run( const std::string& inputFileName, const std::string& maskFileName, const std::string& outputFileName, bool transformToZ, int skipBegin, int skipEnd, float constant )
{
  std::cout << "Read " << inputFileName << std::endl;
  ImageType::Pointer image = LoadImage( inputFileName );
  std::cout << "Read " << maskFileName << std::endl;
  MaskType::Pointer mask = LoadMask( maskFileName );

  ImageType::RegionType region = image->GetLargestPossibleRegion();
  ImageType::SizeType size = region.GetSize();

  int voxels = size[ 0 ] * size[ 1 ] * size[ 2 ];
  int timepoints = size[ 3 ];

  PixelType* data = image->GetPixelContainer()->GetBufferPointer();
  MaskPixelType* dataMask = const_cast< MaskType::ImageType* >( mask->GetImage() )->GetPixelContainer()->GetBufferPointer();

  DataType source( timepoints, voxels, data );
  MatrixType matrix = source.extract( source.rows() - skipBegin - skipEnd, source.cols(), skipBegin, 0 );

  std::vector< int > indices;
  std::vector< int > inverseIndices;
  int rowCount = 0;
  for( int i = 0; i < voxels; ++i )
    {
    if ( dataMask[ i ] )
      {
      indices.push_back( i );
      inverseIndices.push_back( rowCount++ );
      }
    else
      {
      inverseIndices.push_back( -1 );
      }
    }

  int validVoxels = indices.size();
  std::cout << validVoxels << " voxels with " << matrix.rows() << " timepoints" << std::endl;

  OutputImageType::Pointer outputImage = OutputImageType::New();
  OutputImageType::RegionType outputRegion;
  OutputImageType::SizeType outputSize;
  outputSize[ 0 ] = validVoxels;
  outputSize[ 1 ] = validVoxels;
  outputRegion.SetSize( outputSize );
  outputImage->SetRegions( outputRegion );
  outputImage->Allocate();
  outputImage->FillBuffer( itk::NumericTraits< OutputPixelType >::Zero );
  OutputPixelType* outputData = outputImage->GetPixelContainer()->GetBufferPointer();
  OutputDataType output( validVoxels, validVoxels, outputData );

  std::cout << "De-mean" << std::endl;
  std::vector< VectorType > vectors;
  std::vector< PixelType > sums, stds;
  for( int i = 0; i < voxels; ++i )
    {
    if ( !dataMask[ i ] )
      {
      continue;
      }

    VectorType a = matrix.get_column( i );
    a -= a.mean();
    vectors.push_back( a );

    PixelType sum = 0;
    for( int j = a.size() - 1; j >= 0; --j )
      {
      sum += a( j ) * a( j );
      }

    PixelType std = vcl_sqrt( sum / static_cast< PixelType >( a.size() - 1 ) );

    sums.push_back( sum );
    stds.push_back( std );
    }

  image = 0;
  std::cout << "Calculate correlation coefficients" << std::endl;

  const PixelType rowSize = static_cast< PixelType >( matrix.rows() - 1 );
  OutputPixelType factor = 0.5 / vcl_sqrt( 1. / static_cast< OutputPixelType >( matrix.rows() - 3 ) );

  for( int i = 0; i < validVoxels; ++i )
    {
    const VectorType& a = vectors[ i ];

    for( int j = i + 1; j < validVoxels; ++j )
      {
      const VectorType& b = vectors[ j ];
      OutputPixelType r = static_cast< OutputPixelType >( ( dot_product< PixelType >( a, b ) / rowSize ) / ( stds[ i ] * stds[ j ] ) );
      OutputPixelType z = factor * vcl_log( ( 1. + r ) / ( 1. - r ) );

      output( i, j ) = constant + ( transformToZ ? z : r );
      output( j, i ) = constant + ( transformToZ ? z : r );
      }
    }

  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputFileName.c_str() );
  writer->SetInput( outputImage );

  std::cout << "Write " << outputFileName << std::endl;
  writer->Update();
}

RMatrix::ImageType::Pointer RMatrix::LoadImage( const std::string& filename )
{
  typedef itk::ImageFileReader< ImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( filename.c_str() );
  reader->Update();

  return reader->GetOutput();
}

RMatrix::MaskType::Pointer RMatrix::LoadMask( const std::string& filename )
{
  typedef itk::ImageFileReader< MaskImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( filename.c_str() );
  reader->Update();

  MaskType::Pointer mask = MaskType::New();
  mask->SetImage( reader->GetOutput() );
  return mask;
}

} // end namespace fcg


int main( int argc, char ** argv )
{
  tkd::CmdParser p( "rmatrix", "Calculate matrix of voxel-wise correlation coefficients" );

  std::string inputFileName;
  std::string maskFileName;
  std::string outputFileName;

  p.AddArgument( inputFileName, "input" )
    ->AddAlias( "i" )
    ->SetInput( "filename" )
    ->SetDescription( "Input 4D time series image" )
    ->SetRequired( true );

  p.AddArgument( maskFileName, "mask" )
    ->AddAlias( "m" )
    ->SetInput( "filename" )
    ->SetDescription( "Input 3D mask image" )
    ->SetRequired( true );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
    ->SetInput( "filename" )
    ->SetDescription( "Output 2D image: matrix of correlation coefficients" )
    ->SetRequired( true );

  bool transformToZ = false;
  p.AddArgument( transformToZ, "z-transform" )
    ->AddAlias( "z" )
    ->SetInput( "boolean" )
    ->SetDescription( "Output Fisher's z'-transformed r-values" );

  float constant = 0;
  p.AddArgument( constant, "constant" )
      ->AddAlias( "c" )
      ->SetInput( "float" )
      ->SetDescription( "Add constant to all r or z-values" )
	  ->SetRequired( false );

  std::vector< int > skip;
  skip.push_back( 50 );
  skip.push_back( 50 );

  p.AddArgument( skip, "skip" )
    ->AddAlias( "s" )
    ->SetMinMax( 2, 2 )
    ->SetDescription( "Skip samples at begin/end (default: 50 50)" );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  fcg::RMatrix rMatrix;
  rMatrix.Run( inputFileName, maskFileName, outputFileName, transformToZ, skip[ 0 ], skip[ 1 ], constant );

  return 0;
}
