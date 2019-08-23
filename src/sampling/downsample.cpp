#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBSplineDownsampleImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkJoinSeriesImageFilter.h"
#include "tkdCmdParser.h"

int main( int argc, char ** argv )
{
  typedef itk::Image< float, 3 > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typedef itk::Image< float, 2 > SliceType;
  typedef itk::ExtractImageFilter< ImageType, SliceType > ExtractType;
  typedef itk::JoinSeriesImageFilter< SliceType, ImageType > JoinType;
  typedef itk::BSplineDownsampleImageFilter< SliceType, SliceType > FilterType;

  tkd::CmdParser p( "downsample", "Downsample multi-slice dataset slice-by-slice by factor 2" );

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
  ExtractType::Pointer e = ExtractType::New();
  JoinType::Pointer j = JoinType::New();

  r->SetFileName( inputFileName.c_str() );
  r->Update();

  e->SetInput( r->GetOutput() );

  ImageType::RegionType region = r->GetOutput()->GetLargestPossibleRegion();
  ImageType::SizeType size = region.GetSize();
  ImageType::IndexType index = region.GetIndex();

  int slices = size[ 2 ];
  size[ 2 ] = 0;

  region.SetSize( size );

  std::vector< SliceType::Pointer > output;

  for( int i = 0; i < slices; ++i )
    {
    index[ 2 ] = i;
    region.SetIndex( index );

    e = ExtractType::New();
    e->SetInput( r->GetOutput() );
    e->SetExtractionRegion( region );
    e->Update();

    f = FilterType::New();
    f->SetInput( e->GetOutput() );
    f->Update();

    SliceType::Pointer o = f->GetOutput();
    o->DisconnectPipeline();
    f = 0;

    output.push_back( o );

    j->SetInput( i, o );
    }

  j->SetSpacing( r->GetOutput()->GetSpacing()[ 2 ] );

  w->SetFileName( outputFileName.c_str() );
  w->SetInput( j->GetOutput() );
  w->Update();

  return 0;
}
