#include "mrfitCommon.h"
#include "itkMRT1ParameterMap3DImageFilter.h"
#include "itkMRT2ParameterMap3DImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

namespace mrfit
{

struct MRFit::Impl
{
  typedef float PixelType;
  typedef itk::Image< PixelType, 3 > OutputImageType;
  typedef OutputImageType::Pointer OutputImagePointer;

  InputImagePointer Input;
  TimesArray Times;
  std::vector< OutputImagePointer > Output;
};

MRFit::MRFit() : m_Impl( new Impl, true )
{
}

MRFit::~MRFit()
{

}

void MRFit::SetInput( InputImagePointer input )
{
  m_Impl->Input = input;
}

void MRFit::SetTimes( const TimesArray& times )
{
  m_Impl->Times = times;
}

void MRFit::FitT1( T1FittingType type, double maxT1 )
{
	typedef itk::MRT1ParameterMap3DImageFilter< PixelType, PixelType > T1FilterType;
	T1FilterType::Pointer filter = T1FilterType::New();

	T1FilterType::TimeContainerType::Pointer times
		= T1FilterType::TimeContainerType::New();

	InputImagePointer input = m_Impl->Input;
	InputImageType::RegionType inputRegion = input->GetLargestPossibleRegion();
	InputImageType::SizeType inputSize = inputRegion.GetSize();

	typedef T1FilterType::MRImageType ImageType;
	typedef itk::ExtractImageFilter< InputImageType, ImageType > ExtractType;

	unsigned int numberOfTimes = inputSize[ 3 ];
	for( unsigned int i = 0; i < numberOfTimes; ++i ) {
		ExtractType::Pointer extract = ExtractType::New();
		InputImageType::RegionType extractRegion;
		InputImageType::SizeType extractSize;
		InputImageType::IndexType extractIndex;

		extractSize[ 0 ] = inputSize[ 0 ];
		extractSize[ 1 ] = inputSize[ 1 ];
		extractSize[ 2 ] = inputSize[ 2 ];
		extractSize[ 3 ] = 0;

		extractIndex.Fill( 0 );
		extractIndex[ 3 ] = i;

		extractRegion.SetSize( extractSize );
		extractRegion.SetIndex( extractIndex );

		extract->SetInput( input );
		extract->SetExtractionRegion( extractRegion );
		extract->Update();

		ImageType::Pointer time = extract->GetOutput();
		time->DisconnectPipeline();
		extract = 0;

		T1FilterType::TimeType timeType = m_Impl->Times[ i ];

	    filter->AddMRImage( timeType, time );
	}


	filter->SetAlgorithm( type );
	filter->SetMaxT1Time( maxT1 );
	filter->Update();

	m_Impl->Output.clear();

	for( int i = 0; i < 4; ++i )
	{
		typedef itk::VectorIndexSelectionCastImageFilter<
	            itk::VectorImage< PixelType, MapDimension >,
	            itk::Image< PixelType, MapDimension > > VectorIndexSelectionCastImageFilterType;

	    VectorIndexSelectionCastImageFilterType::Pointer extractComponents =
				VectorIndexSelectionCastImageFilterType::New();

	    extractComponents->SetInput( filter->GetOutput() );
	    extractComponents->SetIndex( i );
	    extractComponents->Update();

	    m_Impl->Output.push_back( extractComponents->GetOutput() );
	}
}

void MRFit::FitT2( T2FittingType type, double maxT2 )
{
  typedef itk::MRT2ParameterMap3DImageFilter< PixelType, PixelType > T2FilterType;
  T2FilterType::Pointer filter = T2FilterType::New();

  T2FilterType::EchoTimeContainerType::Pointer times
    = T2FilterType::EchoTimeContainerType::New();

  InputImagePointer input = m_Impl->Input;
  InputImageType::RegionType inputRegion = input->GetLargestPossibleRegion();
  InputImageType::SizeType inputSize = inputRegion.GetSize();

  typedef T2FilterType::MREchoImageType EchoImageType;
  typedef itk::ExtractImageFilter< InputImageType, EchoImageType > ExtractType;

  unsigned int numberOfEchos = inputSize[ 3 ];
  for( unsigned int i = 0; i < numberOfEchos; ++i )
    {
    ExtractType::Pointer extract = ExtractType::New();
    InputImageType::RegionType extractRegion;
    InputImageType::SizeType extractSize;
    InputImageType::IndexType extractIndex;

    extractSize[ 0 ] = inputSize[ 0 ];
    extractSize[ 1 ] = inputSize[ 1 ];
    extractSize[ 2 ] = inputSize[ 2 ];
    extractSize[ 3 ] = 0;

    extractIndex.Fill( 0 );
    extractIndex[ 3 ] = i;

    extractRegion.SetSize( extractSize );
    extractRegion.SetIndex( extractIndex );

    extract->SetInput( input );
    extract->SetExtractionRegion( extractRegion );
    extract->Update();

    EchoImageType::Pointer echo = extract->GetOutput();
    echo->DisconnectPipeline();
    extract = 0;

    T2FilterType::EchoTimeType time = m_Impl->Times[ i ];

    filter->AddMREchoImage( time, echo );
    }

  filter->SetAlgorithm( type );
  filter->SetMaxT2Time( maxT2 );
  //filter->PerformR2MappingOn();
  filter->Update();

  m_Impl->Output.clear();

  for( int i = 0; i < 4; ++i )
    {
    typedef itk::VectorIndexSelectionCastImageFilter<
            itk::VectorImage< PixelType, MapDimension >,
            itk::Image< PixelType, MapDimension > > VectorIndexSelectionCastImageFilterType;

    VectorIndexSelectionCastImageFilterType::Pointer extractComponents =
      VectorIndexSelectionCastImageFilterType::New();

    extractComponents->SetInput( filter->GetOutput() );
    extractComponents->SetIndex( i );
    extractComponents->Update();

    m_Impl->Output.push_back( extractComponents->GetOutput() );
    }
}

MRFit::MapImagePointer MRFit::GetMap( int index )
{
  return m_Impl->Output[ index ];
}

} // end namespace mrfit
