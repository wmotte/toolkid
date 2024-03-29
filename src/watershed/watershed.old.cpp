// The following example illustrates how to preprocess and segment images
// using the \doxygen{WatershedImageFilter}. Note that the care with which
// the data is preprocessed will greatly affect the quality of your result.
// Typically, the best results are obtained by preprocessing the original
// image with an edge-preserving diffusion filter, such as one of the
// anisotropic diffusion filters, or with the bilateral image filter.  As
// noted in Section~\ref{sec:AboutWatersheds}, the height function used as
// input should be created such that higher positive values correspond to
// object boundaries.  A suitable height function for many applications can
// be generated as the gradient magnitude of the image to be segmented.
//
// The \doxygen{VectorGradientMagnitudeAnisotropicDiffusionImageFilter} class
// is used to smooth the image and the
// \doxygen{VectorGradientMagnitudeImageFilter} is used to generate the
// height function.  We begin by including all preprocessing filter header
// files and the header file for the WatershedImageFilter.  We
// use the vector versions of these filters because the input data is a color
// image.
//
//
#include <iostream>

#include "itkVectorGradientAnisotropicDiffusionImageFilter.h"
#include "itkVectorGradientMagnitudeImageFilter.h"
#include "itkWatershedImageFilter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorCastImageFilter.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkScalarToRGBPixelFunctor.h"

int main( int argc, char *argv[] )
{
  if (argc < 8 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage outputImage conductanceTerm diffusionIterations lowerThreshold outputScaleLevel gradientMode " << std::endl;
    return 1;
    }

  // We now declare the image and pixel types to use for instantiation of the
  // filters.  All of these filters expect real-valued pixel types in order to
  // work properly.  The preprocessing stages are done directly on the
  // vector-valued data and the segmentation is done using floating point
  // scalar data.  Images are converted from RGB pixel type to
  // numerical vector type using \doxygen{VectorCastImageFilter}.

  typedef itk::RGBPixel<unsigned char>   RGBPixelType;
  typedef itk::Image<RGBPixelType, 2>    RGBImageType;
  typedef itk::Vector<float, 3>          VectorPixelType;
  typedef itk::Image<VectorPixelType, 2> VectorImageType;
  typedef itk::Image<unsigned long, 2>   LabeledImageType;
  typedef itk::Image<float, 2>           ScalarImageType;

  // The various image processing filters are declared using the types created
  // above and eventually used in the pipeline.

  typedef itk::ImageFileReader<RGBImageType> FileReaderType;
  typedef itk::VectorCastImageFilter<RGBImageType, VectorImageType> CastFilterType;
  typedef itk::VectorGradientAnisotropicDiffusionImageFilter<VectorImageType, VectorImageType>  DiffusionFilterType;
  typedef itk::VectorGradientMagnitudeImageFilter<VectorImageType> GradientMagnitudeFilterType;
  typedef itk::WatershedImageFilter<ScalarImageType> WatershedFilterType;
  typedef itk::ImageFileWriter<RGBImageType> FileWriterType;

  FileReaderType::Pointer reader = FileReaderType::New();
  reader->SetFileName(argv[1]);

  CastFilterType::Pointer caster = CastFilterType::New();

  //
  // Next we instantiate the filters and set their parameters.  The first
  // step in the image processing pipeline is diffusion of the color input
  // image using an anisotropic diffusion filter.  For this class of filters,
  // the CFL condition requires that the time step be no more than 0.25 for
  // two-dimensional images, and no more than 0.125 for three-dimensional
  // images.  The number of iterations and the conductance term will be taken
  // from the command line. See
  // Section~\ref{sec:EdgePreservingSmoothingFilters} for more information on
  // the ITK anisotropic diffusion filters.
  //

  DiffusionFilterType::Pointer diffusion = DiffusionFilterType::New();
  diffusion->SetNumberOfIterations( atoi(argv[4]) );
  diffusion->SetConductanceParameter( atof(argv[3]) );
  diffusion->SetTimeStep(0.125);

  // The ITK gradient magnitude filter for vector-valued images can optionally
  // take several parameters.  Here we allow only enabling or disabling
  // of principal component analysis.
  //

  GradientMagnitudeFilterType::Pointer
    gradient = GradientMagnitudeFilterType::New();

  if ( pca ) {
	  gradient -> SetUsePrincipleComponentsOn();
  } else {
	  gradient -> SetUsePrincipleComponentsOff();
  }


  //
  // Finally we set up the watershed filter.  There are two parameters.
  // \code{Level} controls watershed depth, and \code{Threshold} controls the
  // lower thresholding of the input.  Both parameters are set as a
  // percentage (0.0 - 1.0) of the maximum depth in the input image.
  //

  WatershedFilterType::Pointer watershed = WatershedFilterType::New();
  watershed->SetLevel( atof(argv[6]) );
  watershed->SetThreshold( atof(argv[5]) );

  // The output of WatershedImageFilter is an image of unsigned long integer
  // labels, where a label denotes membership of a pixel in a particular
  // segmented region.  This format is not practical for visualization, so
  // for the purposes of this example, we will convert it to RGB pixels.  RGB
  // images have the advantage that they can be saved as a simple png file
  // and viewed using any standard image viewer software.  The
  // \subdoxygen{Functor}{ScalarToRGBPixelFunctor} class is a special
  // function object designed to hash a scalar value into an
  // \doxygen{RGBPixel}. Plugging this functor into the
  // \doxygen{UnaryFunctorImageFilter} creates an image filter for that
  // converts scalar images to RGB images.

  typedef itk::Functor::ScalarToRGBPixelFunctor<unsigned long>
    ColorMapFunctorType;
  typedef itk::UnaryFunctorImageFilter<LabeledImageType,
    RGBImageType, ColorMapFunctorType> ColorMapFilterType;
  ColorMapFilterType::Pointer colormapper = ColorMapFilterType::New();

  FileWriterType::Pointer writer = FileWriterType::New();
  writer->SetFileName(argv[2]);

  // The filters are connected into a single pipeline, with readers and
  // writers at each end.

  caster->SetInput(reader->GetOutput());
  diffusion->SetInput(caster->GetOutput());
  gradient->SetInput(diffusion->GetOutput());
  watershed->SetInput(gradient->GetOutput());
  colormapper->SetInput(watershed->GetOutput());
  writer->SetInput(colormapper->GetOutput());

  try
    {
    writer->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << e << std::endl;
    }

  return 0;
}

// female head and neck cryosection data.  At left is the original image.  The
// image in the middle was generated with parameters: conductance = 2.0,
// iterations = 10, threshold = 0.0, level = 0.05, principal components = on.
// The image on the right was generated with parameters: conductance = 2.0,
// iterations = 10, threshold = 0.001, level = 0.15, principal components =
// off.
//
// Tuning the filter parameters for any particular application is a process
// of trial and error.  The \emph{threshold} parameter can be used to great
// effect in controlling oversegmentation of the image.  Raising the
// threshold will generally reduce computation time and produce output with
// fewer and larger regions.  The trick in tuning parameters is to consider
// the scale level of the objects that you are trying to segment in the
// image.  The best time/quality trade-off will be achieved when the image is
// smoothed and thresholded to eliminate features just below the desired
// scale.
//
// Figure~\ref{fig:outputWatersheds} shows output from the example code.  The
// input image is taken from the Visible Human female data around the right
// eye.  The images on the right are colorized watershed segmentations with
// parameters set to capture objects such as the optic nerve and
// lateral rectus muscles, which can be seen just above and to the left and
// right of the eyeball.  Note that a critical difference between the two
// segmentations is the mode of the gradient magnitude calculation.
//
// A note on the computational complexity of the watershed algorithm is
// warranted.  Most of the complexity of the ITK implementation lies in
// generating the hierarchy. Processing times for this stage are non-linear
// with respect to the number of catchment basins in the initial segmentation.
// This means that the amount of information contained in an image is more
// significant than the number of pixels in the image.  A very large, but very
// flat input take less time to segment than a very small, but very detailed
// input.

