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
  tkd::CmdParser p( "rotz", "Convert correlation coefficient r-matrix to Fisher's z'-matrix" );

  std::string inputFileName, outputFileName;

  p.AddArgument( inputFileName, "input" )
    ->AddAlias( "i" )
    ->SetInput( "filename" )
    ->SetDescription( "Input image: 2D r-matrix" )
    ->SetRequired( true );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
    ->SetInput( "filename" )
    ->SetDescription( "Output image: 2D z'-matrix" )
    ->SetRequired( true );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  typedef float PixelType;
  typedef itk::Image< PixelType, 2 > ImageType;
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
    it.Set( rs::Common< PixelType >::RtoZ( it.Value() ) );
    }

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputFileName.c_str() );
  writer->SetInput( image );
  writer->Update();

  return 0;
}

