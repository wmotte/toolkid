#ifndef __fidCommon_h__
#define __fidCommon_h__
#include "itkImage.h"
#include "itkPermuteAxesImageFilter.h"
#include "itkFlipImageFilter.h"
#include "vnl/vnl_matrix.h"
#include <vector>
#include "procparser.h"
#include "itkExceptionObject.h"
#include <sstream>

#define fidThrowMacro( text ) \
{ \
  :: std :: ostringstream oss; \
  oss << text; \
  throw :: itk :: ExceptionObject( __FILE__, __LINE__, oss.str().c_str() ); \
}

namespace fid
{

class Common
{
public:
  typedef float PrecisionType;
  typedef std::complex< PrecisionType > ComplexType;
  typedef itk::Image< PrecisionType, 3 > ImageType;
  typedef itk::Image< PrecisionType, 4 > SeriesType;
  typedef itk::Image< ComplexType, 4 > ComplexSeriesType;

  template< class TPixel, unsigned int VDimension >
  static
  typename itk::Image< TPixel, VDimension >::Pointer
  Permute( typename itk::Image< TPixel, VDimension >::Pointer input, int* order )
  {
    typedef itk::Image< TPixel, VDimension > ImageType;
    typedef itk::PermuteAxesImageFilter< ImageType > FilterType;
    typename FilterType::PermuteOrderArrayType permutation;
    for( unsigned int i = 0; i < VDimension; ++i )
      {
      permutation[ i ] = order[ i ];
      }
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput( input );
    filter->SetOrder( permutation );
    filter->Update();

    typename ImageType::Pointer output = filter->GetOutput();
    output->DisconnectPipeline();
    filter = 0;

    typename ImageType::PointType outputOrigin = output->GetOrigin();
    typename ImageType::PointType inputOrigin = input->GetOrigin();

    for( unsigned int i = 0; i < VDimension; ++i )
      {
      outputOrigin[ i ] = inputOrigin[ order[ i ] ];
      }

    output->SetOrigin( outputOrigin );
    return output;
  }

  template< class TPixel, unsigned int InputDimension, unsigned int OutputDimension >
  static
  typename itk::Image< TPixel, OutputDimension >::Pointer
  Reshape( typename itk::Image< TPixel, InputDimension >::Pointer input, int* dimensions, bool manageMemory = false )
  {
    typedef itk::Image< TPixel, InputDimension > InputImageType;
    typedef itk::Image< TPixel, OutputDimension > OutputImageType;

    typename OutputImageType::Pointer output = OutputImageType::New();
    typename OutputImageType::RegionType region;
    typename OutputImageType::SizeType size;

    for( unsigned int i = 0; i < OutputDimension; ++i )
      {
      size[ i ] = dimensions[ i ];
      }

    region.SetSize( size );
    output->SetRegions( region );

    output->GetPixelContainer()->SetImportPointer(
        input->GetPixelContainer()->GetBufferPointer(),
        region.GetNumberOfPixels(),
        manageMemory );

    if ( manageMemory )
      {
      input->GetPixelContainer()->ContainerManageMemoryOff();
      }

    return output;
  }

  template< class TPixel, unsigned int VDimension >
  static
  typename itk::Image< TPixel, VDimension >::Pointer
  Flip( typename itk::Image< TPixel, VDimension >::Pointer input, bool* axes, bool flipAboutOrigin = true )
  {
    typedef itk::Image< TPixel, VDimension > ImageType;
    typedef itk::FlipImageFilter< ImageType > FilterType;
    typename FilterType::FlipAxesArrayType flip;
    for( unsigned int i = 0; i < VDimension; ++i )
      {
      flip[ i ] = axes[ i ];
      }
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput( input );
    filter->SetFlipAxes( flip );
    filter->SetFlipAboutOrigin( flipAboutOrigin );
    filter->Update();

    typename ImageType::Pointer output = filter->GetOutput();
    output->DisconnectPipeline();
    filter = 0;

    output->SetOrigin( input->GetOrigin() );

    return output;
  }

  static ComplexSeriesType::Pointer ReorientCoronal( ComplexSeriesType::Pointer input, const Procparser& pp, bool correctGradientDirection = false );
  static std::vector< PrecisionType > GetSlices( const Procparser& pp );
  static std::vector< int > GetSliceTable( const std::vector< float >& pp );
  static SeriesType::Pointer ComplexToMagnitude( ComplexSeriesType::Pointer input );
  static SeriesType::Pointer ComplexToReal( ComplexSeriesType::Pointer input );
  static SeriesType::Pointer ComplexToImaginary( ComplexSeriesType::Pointer input );
  static std::string GetAbsolutePath( const std::string& execPath, const std::string& fileName );
  static bool LoadTable( const std::string& fileName, vnl_matrix< int >& table );
  static void SaveImage( ImageType::Pointer image, const std::string& filename );
  static void SaveImage( SeriesType::Pointer image, const std::string& filename, int index = 0 );
  static void SaveSeries( SeriesType::Pointer series, const std::string& filename );

};

} // end namespace fid

#endif /*__fidCommon_h__*/
