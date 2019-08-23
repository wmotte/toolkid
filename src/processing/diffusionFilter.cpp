
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"

int main( int argc, char * argv[] )
{
  if( argc < 6 ) 
    { 
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "  input  output ";
    std::cerr << "numberOfIterations (e.g. 5)  timeStep (.e.g 0.0625)  conductance (.e.g. 3 ) useImageSpacingon/off (e.g. 1)" << std::endl;
    return EXIT_FAILURE;
    }
 
  typedef float PixelType; 

  typedef itk::Image< PixelType,  3 >   InputImageType;
  typedef itk::Image< PixelType, 3 >   OutputImageType;
  typedef itk::ImageFileReader< InputImageType >  ReaderType;
  typedef itk::CurvatureAnisotropicDiffusionImageFilter<
               InputImageType, OutputImageType >  FilterType;

  FilterType::Pointer filter = FilterType::New();
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  filter->SetInput( reader->GetOutput() );

  const unsigned int numberOfIterations = atoi( argv[3] );
  const double       timeStep = atof( argv[4] );
  const double       conductance = atof( argv[5] );
  const bool         useImageSpacing = (argc != 6);

  filter->SetNumberOfIterations( numberOfIterations );
  filter->SetTimeStep( timeStep );
  filter->SetConductanceParameter( conductance );
  if (useImageSpacing)
    {
    filter->UseImageSpacingOn();
    }
  filter->Update();

  typedef itk::Image< PixelType, 3 > WriteImageType;

  
  typedef itk::ImageFileWriter< WriteImageType >  WriterType;

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}

