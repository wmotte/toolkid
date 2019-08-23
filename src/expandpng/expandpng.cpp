#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkRGBPixel.h"
#include "tkdCmdParser.h"

int main( int argc, char ** argv )
{
  typedef itk::RGBPixel< unsigned char > PixelType;
  const unsigned int Dimension = 2;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typedef itk::ImageRegionIterator< ImageType > IteratorType;
  typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;
  typedef ImageType::RegionType RegionType;
  typedef ImageType::SizeType SizeType;
  typedef ImageType::IndexType IndexType;

  tkd::CmdParser p( "expandpng", "Combine PNG images" );

  std::vector< std::string > inputFileNames;
  std::string outputFileName;

  p.AddArgument( inputFileNames, "input" )
    ->AddAlias( "i" )
    ->SetInput( "filenames" )
    ->SetDescription( "Input PNG images" )
    ->SetRequired( true );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
    ->SetInput( "filename" )
    ->SetDescription( "Output PNG image" )
    ->SetRequired( true );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  int numberOfImages = inputFileNames.size();
  std::vector< ImageType::Pointer > images;
  std::vector< SizeType > sizes;

  SizeType size;
  size.Fill( 0 );

  for( int i = 0; i < numberOfImages; ++i )
    {
    std::cout << "Read " << inputFileNames[ i ] << std::endl;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( inputFileNames[ i ].c_str() );
    reader->Update();

    ImageType::Pointer image = reader->GetOutput();
    image->DisconnectPipeline();
    images.push_back( image );

    RegionType region = image->GetLargestPossibleRegion();
    sizes.push_back( region.GetSize() );

    size[ 0 ] = size[ 0 ] < sizes[ i - 2 ][ 0 ] ? sizes[ i - 2 ][ 0 ] : size[ 0 ];
    size[ 1 ] += sizes[ i - 2 ][ 1 ];
    }

  std::cout << "Output: " << size << std::endl;

  ImageType::Pointer output = ImageType::New();
  output->CopyInformation( images[ 0 ] );

  RegionType region;
  region.SetSize( size );

  output->SetRegions( region );
  output->Allocate();

  IndexType index;
  index.Fill( 0 );


  for( unsigned int i = 0; i < images.size(); ++i )
    {
    region.SetSize( sizes[ i ] );
    region.SetIndex( index );

    ConstIteratorType inputIt( images[ i ], images[ i ]->GetLargestPossibleRegion() );
    for( IteratorType it( output, region ); !it.IsAtEnd(); ++it, ++inputIt )
      {
      it.Set( inputIt.Get() );
      }

    index[ 1 ] += sizes[ i ][ 1 ];
    }

  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( output );
  writer->SetFileName( outputFileName.c_str() );
  writer->Update();

  return 0;
}
