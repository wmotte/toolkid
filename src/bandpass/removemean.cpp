#include "itkImage.h"
#include "vnl/vnl_matrix_ref.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "vnl/vnl_vector.h"
#include "tkdCmdParser.h"

int main( int argc, char ** argv )
{
  tkd::CmdParser p( "removemean", "Remove mean signal fluctuations from 4-D dataset" );

  std::string inputFileName, outputFileName, maskFileName;

  p.AddArgument( inputFileName, "input" )
    ->AddAlias( "i" )
    ->SetInput( "filename" )
    ->SetDescription( "4D input image" )
    ->SetRequired( true );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
    ->SetInput( "filename" )
    ->SetDescription( "4D filtered output image" )
    ->SetRequired( true );

  p.AddArgument( maskFileName, "mask" )
    ->AddAlias( "m" )
    ->SetDescription( "3-D mask" );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }


  typedef itk::Image< float, 4 > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typedef itk::Image< unsigned char, 3 > MaskType;
  typedef itk::ImageFileReader< MaskType > MaskReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  MaskReaderType::Pointer maskReader = MaskReaderType::New();

  reader->SetFileName( inputFileName.c_str() );
  reader->Update();

  ImageType::Pointer image = reader->GetOutput();
  reader = 0;

  float* buffer = image->GetPixelContainer()->GetBufferPointer();
  int timepoints = image->GetLargestPossibleRegion().GetSize()[ 3 ];
  int voxels = image->GetLargestPossibleRegion().GetNumberOfPixels() / timepoints;

  vnl_matrix_ref< float > matrix( timepoints, voxels, buffer );

  vnl_vector< int > mask( voxels );
  mask.fill( 1 );

  if ( maskFileName != "" )
    {
    maskReader->SetFileName( maskFileName.c_str() );
    maskReader->Update();
    unsigned char* maskBuffer = maskReader->GetOutput()->GetPixelContainer()->GetBufferPointer();

    for( int i = 0; i < voxels; ++i )
      {
      mask[ i ] = ( maskBuffer[ i ] != 0 ) ? 1 : 0;
      }
    }

  vnl_vector< float > mean( timepoints );
  mean.fill( 0 );
  int count = 0;
  for( int i = 0; i < voxels; ++i )
    {
    if ( mask[ i ] == 0 )
      {
      continue;
      }

    ++count;
    mean += matrix.get_column( i );
    }

  mean /= static_cast< float >( count );

  for( int i = 0; i < voxels; ++i )
    {
    matrix.set_column( i, matrix.get_column( i ) - mean );
    }

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputFileName.c_str() );
  writer->SetInput( image );
  writer->Update();

  return 0;
}

