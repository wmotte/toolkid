#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRGBPixel.h"
#include "itkScalarToRGBColormapImageFilter.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageSliceConstIteratorWithIndex.h"
#include "itkExtractImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkCustomColormapFunctor.h"
#include <fstream>
#include "itkShiftScaleImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkPermuteAxesImageFilter.h"
#include "itkFlipImageFilter.h"

#include "graphCommon.h"

namespace colormap
{

/**
 * Option list.
 */
struct parameters
{
  std::vector< std::string > inputFileNames;
  std::string backgroundFileName;
  std::string outputFileName;
  std::string maskFileName;
  int sliceDirection;
  bool noMask;
  int interpolationOrder;
  double percentile;
  double scaleFactor;

  std::vector< double > overlaysLowerThreshold;
  std::vector< double > overlaysUpperThreshold;
  std::vector< double > overlayTransparencies;

  double backgroundTransparency;

  double backgroundMin;
  double backgroundMax;

  std::vector< std::string > LUTS;

  int numberOfColumns;

  std::vector< int > selectedSlices;
  int startSlice;
  int endSlice;
  int skipSlices;

  bool addColorBars;
  bool colorBarHorizontal;
  int colorBarWidth;
  
  int defaultR;
  int defaultB;
  int defaultG;
};

/**
 * ColorMap implementation.
 */
class ColorMap
{
public:
  typedef double PixelType;
  typedef itk::RGBPixel< unsigned char > RGBPixelType;

  typedef std::vector< PixelType > VectorType;

  typedef itk::Image< PixelType, 3 > ImageType;
  typedef itk::Image< RGBPixelType, 3 > RGBImageType;

  typedef itk::MinimumMaximumImageFilter< ImageType > MinimumMaximumImageFilterType;
  typedef itk::RescaleIntensityImageFilter< ImageType > RescaleIntensityImageFilterType;
  typedef itk::ThresholdImageFilter< ImageType > ThresholdImageFilterType;
  typedef itk::ImageFileReader< ImageType > ImageFileReaderType;

  typedef itk::ScalarToRGBColormapImageFilter< ImageType, RGBImageType > RGBFilterType;
  typedef itk::Functor::CustomColormapFunctor< PixelType, RGBPixelType > ColormapType;

  typedef itk::Image< RGBPixelType, 2 > OverlayImageType;
  typedef itk::ImageFileWriter< OverlayImageType > OverlayImageFileWriterType;

  typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
  typedef itk::ImageRegionIteratorWithIndex< RGBImageType > RGBIteratorType;
  typedef itk::ImageSliceConstIteratorWithIndex< RGBImageType > ImageSliceConstIteratorWithIndexType;
  typedef itk::FlipImageFilter< OverlayImageType > FlipImageFilterType;
  typedef itk::ExtractImageFilter< ImageType, ImageType > ExtractImageFilterType;
  typedef itk::MaskImageFilter< RGBImageType, ImageType, RGBImageType > MaskImageFilterType;
  typedef itk::IntensityWindowingImageFilter< ImageType, ImageType > IntensityWindowingImageFilterType;

  typedef itk::ImageToImageFilter< ImageType, RGBImageType > ImageToRGBImageFilterType;

  typedef itk::ResampleImageFilter< ImageType, ImageType > FilterType;
  typedef itk::BSplineInterpolateImageFunction< ImageType, double > InterpolatorType;
  typedef itk::PermuteAxesImageFilter< RGBImageType > PermuteFilterType;
  typedef itk::FlipImageFilter< RGBImageType > FlipFilterType;

private:

  std::vector< ImageType::Pointer > m_Images;
  std::vector< RGBImageType::Pointer > m_RGBImages;
  std::vector< RGBImageType::Pointer > m_LUTImages;
  bool m_UseMask;
  bool m_UseBackground;

  ImageType::Pointer m_Mask;
  ImageType::Pointer m_Background;

public:
  /**
   * Run.
   */
  void Run( const parameters& list );

protected:
  /**
   * Permute merged RGB image such that slices are in the Z-direction
   */
  RGBImageType::Pointer PermuteImage( RGBImageType::Pointer image, unsigned int sliceDirection );
  std::vector< int > CalculateSlices( RGBImageType::Pointer image, int skip = 1, int start = -1, int end = -1 );

  /**
   * Rescale background image
   */
  void LevelSetBackground( PixelType min, PixelType max, double percentile = 99.5 );

  /**
   * Merge overlays with background image.
   */
  RGBImageType::Pointer MergeBackground( PixelType transparency, int defaultR = 0, int defaultG = 0, int defaultB = 0 );
  void Percentile( ImageType::Pointer input, ImageType::Pointer mask, double percentile, double& lower, double& upper );
  ImageType::Pointer MaskImage( ImageType::Pointer input, ImageType::Pointer mask );

  /**
   * Crop all images using mask (if specified).
   */
  void CropAndMask( bool applyMask = true );

  /**
   * Flip along axis.
   */
  OverlayImageType::Pointer Flip( const OverlayImageType::Pointer& overlay, bool x = false, bool y = false );

  /**
   * New overlay.
   */
  OverlayImageType::Pointer Overlay(
      unsigned int nCols,
      std::vector< int > selectedSlices,
      const RGBImageType::Pointer& image3D
      );

  /**
   * Write overlay image type to file.
   */
  void WriteOverlay( const OverlayImageType::Pointer& image, const std::string& outputFileName );

  /**
   * Merge all RGB images into one.
   */
  RGBImageType::Pointer MergeRGBImages( const RGBImageType::Pointer& background,
      const VectorType& overlayTransparencies );

  /**
   * Colormap all input images according specified LUTS.
   */
  void GenerateRGB( const std::vector< std::string >& LUTS );
  OverlayImageType::Pointer AddColorBars( OverlayImageType::Pointer input, const std::vector< std::string >& LUTS, const VectorType& lowerThreshold, const VectorType& upperThreshold, int width, bool horizontal );

  /**
   * Return custom colormap if .txt file is specified in application directory.
   */
  ColormapType::Pointer GetCustomColorMap( const std::string& name );

  /**
   * Generate RGB image from scalar input image.
   */
  RGBImageType::Pointer GenerateRGB( const std::string& colormapString, ImageType::Pointer image );

  /**
   * Specify LUTS for all images.
   * First image is considered being background.
   */
  std::vector< std::string > SelectLuts( unsigned int total );

  /**
   * Rescale all images between 0 and 255.
   */
  void RescaleImages( VectorType lower, VectorType upper );

  ImageType::Pointer ResampleImage( ImageType::Pointer image, ImageType::Pointer reference, int interpolationOrder = 3 );

  void ResampleImages( ImageType::Pointer reference, unsigned int sliceDirection, double scaleFactor = 1.0, int interpolationOrder = 3 );

  /**
   * Read input images.
   */
  void ReadImages( const std::vector< std::string >& inputFileNames );

  /**
   * Read mask.
   */
  void ReadMask( const std::string& maskFileName );

  /**
   * Read background.
   */
  void ReadBackground( const std::string& backgroundFileName );

  /**
   * Get minimal voxel intensities.
   */
  VectorType GetMinValues();

  /**
   * Get maximal voxel intensities.
   */
  VectorType GetMaxValues();

  /**
   * Check parameters, and exist if not correct.
   */
  void Checks( const parameters& list );
};

} // end namespace colormap
