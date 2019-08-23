#ifndef __looklocker_h__
#define __looklockerm_h__


#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkPasteImageFilter.h"
#include "itkExtractImageFilter.h"


#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"

#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include "boost/type.hpp"

#include "fidFIDReader.h"
#include "fidFID.h"
#include "fidFourier.h"
#include "fidCommon.h"

#include "tkdCmdParser.h"

#include "procparser.h"

#include "mrfitCommon.h"

namespace mapping
{

	typedef fid::FID::PrecisionType PixelType;
	typedef fid::FID::Pointer FIDPointerType;

	typedef itk::Image< PixelType, 3 > VolumeType;
	typedef itk::Image< PixelType, 4 > TimeSeriesType;
	typedef itk::ImageFileReader< TimeSeriesType > ReaderType;
	typedef itk::ImageFileWriter< VolumeType > WriterType;
	typedef itk::RegionOfInterestImageFilter< TimeSeriesType, TimeSeriesType > RegionOfInterestImageFilterType;
	typedef itk::PasteImageFilter< VolumeType, VolumeType, VolumeType > PasteImageFilterType;
	typedef itk::ExtractImageFilter< TimeSeriesType, VolumeType > ExtractVolumeFilterType;

	typedef std::vector< double > VectorType;
	typedef std::vector< int > SliceTableType;

	/**
	 * Command line argument structure.
	 */
	struct cmdParameters
	{
		std::string inputFileName;
		std::string fidFileName;
		std::string outputFileName;
		std::string outputAFileName;
		std::string outputBFileName;
		float maxT1;
	};

	/**
	 * T1 mapping look locker.
	 */
	class T1MappingLookLocker
	{

	public:

		/**
		 * Run.
		 */
		void Run( const cmdParameters& args );

		/**
		 * Extra function to use from procfid applications.
		 */
		void Process( const TimeSeriesType::Pointer& input, const FIDPointerType& fid, const std::string& outputFileName, const std::string& outputAFileName, const std::string& outputBFileName, PixelType maxT1 );

	private:

		/**
		 * Write output to file.
		 */
		void Write( const VolumeType::Pointer& input, const std::string& fileName );

		/**
		 * Return list of TR's for given slice.
		 */
		VectorType GetTRs( unsigned int slice, unsigned int z, unsigned int t, PixelType trimage, PixelType ti );

		/**
		 * Allocate 3D input, using the 4D input image information.
		 */
		VolumeType::Pointer AllocateVolume( const TimeSeriesType::Pointer& input );

		/**
		 * Return 4D input data.
		 */
		TimeSeriesType::Pointer GetTimeSeries( const std::string& inputFileName );

		/**
		 * Return all time-slices from one given slice.
		 */
		TimeSeriesType::Pointer GetSliceTimeSeries( const unsigned int slice, TimeSeriesType::Pointer input );

		/**
		 * Return FID pointer.
		 */
		FIDPointerType GetFID( const std::string& inputFileName );

		/**
		 * Add slice to map.
		 */
		void AddSliceToVolume( const VolumeType::Pointer& input, VolumeType::Pointer& output, unsigned int slice );

		/**
		 * Return number of Slices.
		 */
		unsigned int GetNumberOfSlices( const TimeSeriesType::Pointer& input );

		/**
		 * Return number of Volumes.
		 */
		unsigned int GetNumberOfVolumes( const TimeSeriesType::Pointer& input );

	}; // class T1MappingLookLocker

} // mapping namespace

#endif // __looklocker_h__
