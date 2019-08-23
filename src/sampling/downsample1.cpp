#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBSplineDownsampleImageFilter.h"

int main( int argc, char ** argv )
{
  typedef itk::Image< float, 3 > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typedef itk::BSplineDownsampleImageFilter< ImageType, ImageType > FilterType;
  
  ReaderType::Pointer r = ReaderType::New();
  WriterType::Pointer w = WriterType::New();
  FilterType::Pointer f = FilterType::New();
  
  r->SetFileName( argv[ 1 ] );
  w->SetFileName( argv[ 2 ] );
  f->SetInput( r->GetOutput() );
  w->SetInput( f->GetOutput() );
  
  w->Update();
  
  return 0;
}
