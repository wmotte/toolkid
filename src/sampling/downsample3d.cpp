#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBSplineDownsampleImageFilter.h"
#include "tkdCmdParser.h"

int main( int argc, char ** argv )
{
  typedef itk::Image< float, 3 > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typedef itk::BSplineDownsampleImageFilter< ImageType, ImageType > FilterType;

  tkd::CmdParser p( "downsample", "Downsample 3D dataset in each direction by a factor 2" );

  std::string inputFileName, outputFileName;

  p.AddArgument( inputFileName, "input" )
    ->AddAlias( "i" )
    ->SetInput( "filename" )
    ->SetDescription( "Input filename" )
    ->SetRequired( true );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
    ->SetInput( "filename" )
    ->SetDescription( "Output filename" )
    ->SetRequired( true );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  ReaderType::Pointer r = ReaderType::New();
  WriterType::Pointer w = WriterType::New();
  FilterType::Pointer f = FilterType::New();

  r->SetFileName( inputFileName.c_str() );
  r->Update();

  f->SetInput( r->GetOutput() );

  w->SetFileName( outputFileName.c_str() );
  w->SetInput( f->GetOutput() );
  w->Update();

  return 0;
}
