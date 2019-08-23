#include "itkImage.h"
#include "vnl/vnl_matrix_ref.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "vnl/vnl_math.h"
#include "tkdCmdParser.h"
#include "rsCommon.h"

int main( int argc, char ** argv )
{
  tkd::CmdParser p( "rmaptoz", "Convert between correlation coefficient r and Fisher z'" );

  std::string inputFileName, outputFileName;
  bool zToR = false;

  p.AddArgument( inputFileName, "input" )
    ->AddAlias( "i" )
    ->SetDescription( "Input image: 4D z' (or r)-map" )
    ->SetRequired( true );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
    ->SetDescription( "Output image: 4D r (or z')-map" )
    ->SetRequired( true );

  p.AddArgument( zToR, "z-to-r" )
    ->AddAlias( "z" )
    ->SetDescription( "Convert input z'-map back to r-values" );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  typedef float PixelType;
  typedef itk::Image< PixelType, 4 > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typedef vnl_matrix_ref< PixelType > MatrixType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputFileName.c_str() );
  reader->Update();

  ImageType::Pointer image = reader->GetOutput();
  reader = 0;

  itk::ImageRegionIterator< ImageType > it( image, image->GetLargestPossibleRegion() );
  for( ; !it.IsAtEnd(); ++it )
    {
    it.Set( zToR ? rs::Common< PixelType >::ZtoR( it.Value() )
                 : rs::Common< PixelType >::RtoZ( it.Value() ) );
    }

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputFileName.c_str() );
  writer->SetInput( image );
  writer->Update();

  return 0;
}

