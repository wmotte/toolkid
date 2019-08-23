#include "fidFIDReader.h"
#include "fidFID.h"
#include "fidFourier.h"
#include "fidCommon.h"
#include "tkdCmdParser.h"
#include "mrfitCommon.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkPasteImageFilter.h"

#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include "boost/type.hpp"

#include "vnl/vnl_matrix.h"

#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
//-----------------------------------------------------------------------------------------

typedef fid::FID::PrecisionType PrecisionType;
typedef fid::FID::ComplexType ComplexType;
typedef fid::FID::DataType DataType;
typedef fid::Fourier< PrecisionType > FourierType;

typedef itk::Image< float, 2 > SliceType;
typedef itk::Image< float, 3 > OutputMapType;
typedef itk::Image< float, 4 > TimeSeriesType;

typedef itk::Image< ComplexType, 4 > OutputImageType;
typedef itk::ImageFileWriter< TimeSeriesType > TimeSeriesWriter;
typedef itk::PermuteAxesImageFilter< TimeSeriesType > PermuteFilter;
typedef itk::RegionOfInterestImageFilter< TimeSeriesType, TimeSeriesType > RegionOfInterestImageFilterType;
typedef itk::PasteImageFilter< OutputMapType, OutputMapType, OutputMapType > PasteImageFilterType;

typedef vnl_matrix< ComplexType > MatrixType;

typedef std::vector< double > VectorType;

//-----------------------------------------------------------------------------------------

std::vector< int > zeroFill;

/****************************************************************************************/

/**
 * Get empty maps to fill with slices...
 */
void GetEmptyMaps( const unsigned int zeroFillPE, const unsigned int zeroFillRO, const std::vector< float >& pss,
		const std::vector< int >& slices, const fid::FID::Pointer fid, OutputMapType::Pointer& t1map, OutputMapType::Pointer& amap,
		OutputMapType::Pointer& bmap, OutputMapType::Pointer& r2map ) {

	t1map = OutputMapType::New();
	amap = OutputMapType::New();
	bmap = OutputMapType::New();
	r2map = OutputMapType::New();

	unsigned int ro = fid -> GetNumberOfPoints();
	unsigned int pe = fid -> Procpar< int > ( "nv" );
	unsigned int ns = fid -> Procpar< int > ( "ns" );

	OutputMapType::RegionType region;
	OutputMapType::SizeType size;
	OutputMapType::IndexType index;

	OutputMapType::PointType origin;
	OutputMapType::SpacingType spacing;

	if ( zeroFill.size() == 2 ) {
		ro = zeroFill[0];
		pe = zeroFill[1];
	}

	index[0] = 0;
	index[1] = 0;
	index[2] = 0;

	// pe x slices x ro
	size[0] = pe;
	size[1] = ns;
	size[2] = ro;

	region.SetIndex( index );
	region.SetSize( size );

	// begin spacing
	float lro = fid -> Procpar< float > ( "lro" );
	float lpe = fid -> Procpar< float > ( "lpe" );
	float sum = 0;
	for ( unsigned int i = 1; i < pss.size(); ++i ) {
		sum += ( pss[slices[i]] - pss[slices[i - 1]] );
	}
	float slcdist = pss.size() == 1 ? .1 : ( sum / static_cast< float > ( pss.size() - 1 ) );
	spacing[0] = lpe / static_cast< float > ( zeroFillPE ) * 10.0;
	spacing[1] = slcdist * 10.0;
	spacing[2] = lro / static_cast< float > ( zeroFillRO ) * 10.0;
	// end spacing

	// begin origin
	float pss0 = fid -> Procpar< float > ( "pss0" ) * 10.0;
	for ( int i = 0; i < 3; ++i ) {
		origin[i] = -0.5 * spacing[i] * static_cast< float > ( size[i] );
	}
	origin[1] += pss0;
	// end origin

	t1map -> SetRegions( region );
	t1map -> SetSpacing( spacing );
	t1map -> SetOrigin( origin );

	amap -> CopyInformation( t1map );
	bmap -> CopyInformation( t1map );
	r2map -> CopyInformation( t1map );

	t1map -> Allocate();
	amap -> Allocate();
	bmap -> Allocate();
	r2map -> Allocate();
}

void getOutputImage( const unsigned int x, const unsigned int y, const unsigned int z, const unsigned int t, const unsigned int zeroFillPE,
		const unsigned int zeroFillRO, const std::vector< float >& pss, const std::vector< int >& slices, const fid::FID::Pointer& fid,
		OutputImageType::Pointer& outputImage, OutputImageType::SizeType& size, OutputImageType::IndexType& index,
		OutputImageType::RegionType& region ) {
	outputImage = OutputImageType::New();

	size[0] = x;
	size[1] = y;
	size[2] = z;
	size[3] = t;

	index[0] = 0;
	index[1] = 0;
	index[2] = 0;
	index[3] = 0;

	region.SetSize( size );
	region.SetIndex( index );

	outputImage -> SetRegions( region );
	outputImage -> Allocate();

	// calculate spacing and field of view

	float lro = fid -> Procpar< float > ( "lro" );
	float lpe = fid -> Procpar< float > ( "lpe" );

	float sum = 0;

	for ( unsigned int i = 1; i < pss.size(); ++i ) {
		sum += ( pss[slices[i]] - pss[slices[i - 1]] );
	}

	float slcdist = pss.size() == 1 ? .1 : ( sum / static_cast< float > ( pss.size() - 1 ) );

	OutputImageType::SpacingType spacing;
	spacing[0] = lpe / static_cast< float > ( zeroFillPE ) * 10.0;
	spacing[1] = lro / static_cast< float > ( zeroFillRO ) * 10.0;
	spacing[2] = slcdist * 10.0;
	spacing[3] = 1;

	outputImage -> SetSpacing( spacing );

	float pss0 = fid -> Procpar< float > ( "pss0" ) * 10.0;

	OutputImageType::PointType origin;
	for ( int i = 0; i < 3; ++i ) {
		origin[i] = -0.5 * spacing[i] * static_cast< float > ( size[i] );
	}
	origin[2] += pss0;
	origin[3] = 0;

	outputImage -> SetOrigin( origin );
}

/**
 * Write 4D-series to disk.
 */
void writeSeriesToDisk( TimeSeriesType::Pointer input, const std::string& outputName ) {

	TimeSeriesWriter::Pointer writer = TimeSeriesWriter::New();
	writer -> SetInput( input );
	writer -> SetFileName( outputName.c_str() );

	try {
		writer -> Update();
	} catch ( itk::ExceptionObject& e ) {
		std::cerr << "Error writing to: " << outputName << std::endl;
		std::cerr << e.GetDescription() << std::endl;
	}
}

/**
 * Return all time-slices from one given slice.
 */
TimeSeriesType::Pointer getTimeSeries( const unsigned int slice, TimeSeriesType::Pointer input ) {

	TimeSeriesType::RegionType inputRegion = input -> GetLargestPossibleRegion();
	TimeSeriesType::SizeType inputSize = inputRegion.GetSize();
	TimeSeriesType::IndexType inputIndex = inputRegion.GetIndex();

	TimeSeriesType::IndexType index;
	TimeSeriesType::SizeType size;
	TimeSeriesType::RegionType region;

	index[0] = inputIndex[0];
	index[1] = slice;
	index[2] = inputIndex[2];
	index[3] = inputIndex[3];

	size[0] = inputSize[0];
	size[1] = 1;
	size[2] = inputSize[2];
	size[3] = inputSize[3];

	region.SetIndex( index );
	region.SetSize( size );

	RegionOfInterestImageFilterType::Pointer filter = RegionOfInterestImageFilterType::New();
	filter -> SetRegionOfInterest( region );
	filter -> SetInput( input );
	filter -> Update();

	TimeSeriesType::Pointer output = filter -> GetOutput();

	return output;
}

/**
 * Construct TR-vector for a specific slice, using BOOST reshape.
 */
VectorType buildSliceTimeTable( const unsigned int slice, const unsigned int ns, const unsigned int ne, const double trimage,
		const double ti ) {

	unsigned int total = ns * ne;

	VectorType vector;
	double queue[total];

	// First get list of all possibilities (600)...
	for ( unsigned int i = 0; i < total; i++ ) {

		double value = ti + i * trimage;
		queue[i] = value;
	}

	typedef boost::multi_array< double, 2 > array;
	typedef boost::const_multi_array_ref< double, 2 > const_array_ref;
	typedef boost::general_storage_order< 2 > storage;

	// ordering similar to matlab's reshape...
	array::size_type ordering[] = { 0, 1 };
	bool ascending[] = { true, true };

	// reshape
	boost::array< array::size_type, 2 > dims = { { ns, ne } };
	const_array_ref matrix( queue, dims, storage( ordering, ascending ) );
	matrix.reshape( dims );

	for ( unsigned int i = 0; i < ne; i++ ) {

		vector.push_back( matrix[slice][i] );

	}

	return vector;
}

/*
 * Add slice to map. Input expects: PE x slices x RO!! (e.g. 64 x 25 x 64).
 */
void addSliceToMaps( const unsigned int slice, mrfit::MRFit::Pointer& fit, OutputMapType::Pointer& t1map, OutputMapType::Pointer& amap,
		OutputMapType::Pointer& bmap, OutputMapType::Pointer& r2map ) {
	// index...
	PasteImageFilterType::InputImageIndexType index;
	index[0] = 0;
	index[1] = slice;
	index[2] = 0;

	// region...
	OutputMapType::RegionType region = fit -> GetMap( mrfit::MRFit::T1_MAP ) -> GetLargestPossibleRegion();

	// filter [ T1 ]
	PasteImageFilterType::Pointer filter1 = PasteImageFilterType::New();
	filter1 -> SetSourceRegion( region );
	filter1 -> SetDestinationIndex( index );

	filter1 -> SetSourceImage( fit -> GetMap( mrfit::MRFit::T1_MAP ) );
	filter1 -> SetDestinationImage( t1map );
	filter1 -> Update();
	t1map = filter1 -> GetOutput();
}

/**
 * Construct T1 map.
 */
void calculateT1( const std::string& outputMapName, const itk::Image< float, 4 >::Pointer input, const fid::FID::Pointer fid,
		const double maxT1, const int algorithm, const std::string& outputA, const std::string& outputB, const std::string& outputR2,
		const int zeroFillPE, const int zeroFillRO, const std::vector< float >& pss, const std::vector< int >& slices ) {
	unsigned int ns = fid -> Procpar< int > ( "ns" );
	unsigned int ne = fid -> Procpar< int > ( "ne" );
	double ti = fid -> Procpar< double > ( "ti" );
	double trimage = fid -> Procpar< double > ( "trimage" );
	std::vector< int > sliceTable = fid::Common::GetSliceTable( fid::Common::GetSlices( fid->GetProcpar() ) );

	OutputMapType::Pointer t1map;
	OutputMapType::Pointer amap;
	OutputMapType::Pointer bmap;
	OutputMapType::Pointer r2map;

	getEmptyMaps( zeroFillPE, zeroFillRO, pss, slices, fid, t1map, amap, bmap, r2map );

	// for each slice...
	for ( unsigned int z = 0; z < ns; z++ ) {

		// get TRs...
		VectorType table = buildSliceTimeTable( sliceTable[z], ns, ne, trimage, ti );

		// get slice time-series...
		TimeSeriesType::Pointer timeSeries = getTimeSeries( z, input );

		// fit slice...
		mrfit::MRFit::Pointer fit = mrfit::MRFit::New();
		fit -> SetInput( timeSeries ); // in seconds...
		fit -> SetTimes( table );
		fit -> FitT1( static_cast< mrfit::MRFit::T1FittingType > ( algorithm ), maxT1 );

		// insert fitted map into empty output maps...
		addSliceToMaps( z, fit, t1map, amap, bmap, r2map );
	}

	// output
	if ( outputMapName != "" ) {
		fid::Common::SaveImage( t1map, outputMapName );
	}
}

/**
 * Set own priority.
 */
void setPriority( const int priority ) {
	setpriority( PRIO_PROCESS, getpid(), priority );
	//std::cout << "Resetted process priority to: " << getpriority( PRIO_PROCESS, getpid() ) << std::endl;
}

/**
 * Reset data to coronal slices (do not use Kajo's ReorientCoronal as this flips the x direction...
 */
OutputImageType::Pointer reorientCoronal( OutputImageType::Pointer input ) {

	//bool flip[] = { true, false, false, false };
	bool flip[] = { false, false, false, false };
	int permutation[] = { 0, 2, 1, 3 };

	input = fid::Common::Permute< ComplexType, 4 >( fid::Common::Flip< ComplexType, 4 >( input, flip, false ), permutation );

	OutputImageType::DirectionType direction;
	direction.SetIdentity();

	input->SetDirection( direction );

	return input;

}

/***************************************************************************************
 * Main.
 ****************************************************************************************/
int main( int argc, char ** argv ) {
	// set priority low...
	setPriority( 15 );

	tkd::CmdParser p( "proclooklocker", "Fourier transform Varian Look-Locker FID" );

	std::string inputFileName;
	std::string outputFileName;
	std::string outputRealFileName;
	std::string outputImaginaryFileName;

	std::string outputA;
	std::string outputB;
	std::string outputR2;

	std::string outputFileNameT1;
	int algorithm = mrfit::MRFit::ABSOLUTE_INVERSION_RECOVERY;

	double maxT1 = 10;
	bool outputKSpace = false;

	p.AddArgument( inputFileName, "input" ) ->AddAlias( "i" ) ->SetDescription( "FID input folder" ) ->SetRequired( true );

	p.AddArgument( outputFileName, "output" ) ->AddAlias( "o" ) ->SetDescription( "Output image" ) ->SetRequired( true );

	p.AddArgument( outputRealFileName, "output-real" ) ->AddAlias( "or" ) ->SetDescription( "Output real image" );

	p.AddArgument( outputImaginaryFileName, "output-imag" ) ->AddAlias( "oi" ) ->SetDescription( "Output imaginary image" );

	p.AddArgument( outputFileNameT1, "output-T1map" ) ->AddAlias( "t1map" ) ->SetDescription( "Perform T1-mapping" );

	p.AddArgument( zeroFill, "zerofill" ) ->AddAlias( "zf" ) ->SetInput( "<int int>" ) ->SetDescription( "Zero-fill matrix" ) ->SetMinMax(
			2, 2 );

	p.AddArgument( maxT1, "maximum-T1" ) ->AddAlias( "max" ) ->SetDescription( "Maximum T1 in seconds (default: 10)" );

	p.AddArgument( algorithm, "algorithm" ) ->AddAlias( "a" ) ->SetDescription( "Fitting algorithm: "
		"0=IDEAL_STEADY_STATE, "
		"1=INVERSION_RECOVERY, "
		"2=ABSOLUTE_INVERSION_RECOVERY, "
		"3=LOOK_LOCKER, "
		"4=ABSOLUTE_LOOK_LOCKER, "
		"5=HYBRID_STEADY_STATE_3PARAM, "
		"6=INVERSION_RECOVERY_3PARAM, "
		"7=ABSOLUTE_INVERSION_RECOVERY_3PARAM "
		"(default: 2)" );

	if ( !p.Parse( argc, argv ) ) {
		p.PrintUsage( std::cout );
		return -1;
	}

	typedef fid::FID::PrecisionType PrecisionType;
	typedef fid::FID::ComplexType ComplexType;
	typedef fid::FID::DataType DataType;
	typedef fid::Fourier< PrecisionType > FourierType;
	typedef vnl_matrix< ComplexType > MatrixType;
	typedef itk::Image< ComplexType, 4 > OutputImageType;

	// Read FID and procpar

	fid::FIDReader::Pointer reader = fid::FIDReader::New();
	reader -> SetFileName( inputFileName );

	try {
		reader -> Read();
	} catch ( itk::ExceptionObject& e ) {
		std::cerr << "Error reading FID: " << e.GetDescription() << std::endl;
		return -1;
	}

	fid::FID::Pointer fid = reader->GetFID();

	// determine dimensions

	int numberOfRO = fid->GetNumberOfPoints();
	int numberOfPE = fid->Procpar< int > ( "nv" );
	int numberOfEchos = fid->Procpar< int > ( "ne" );
	int numberOfZ = fid->Procpar< int > ( "ns" ); // slices


	// reshape

	typedef itk::Image< ComplexType, 4 > ShapedDataType;
	ShapedDataType::Pointer shaped = ShapedDataType::New();
	ShapedDataType::RegionType shapedRegion;
	ShapedDataType::IndexType shapedIndex;
	ShapedDataType::SizeType shapedSize;

	shapedSize[0] = numberOfRO;
	shapedSize[1] = numberOfZ;
	shapedSize[2] = numberOfEchos;
	shapedSize[3] = numberOfPE;

	shapedIndex.Fill( 0 );

	shapedRegion.SetIndex( shapedIndex );
	shapedRegion.SetSize( shapedSize );

	shaped-> SetRegions( shapedRegion );
	shaped-> GetPixelContainer() -> SetImportPointer( fid -> GetData() -> GetPixelContainer() -> GetBufferPointer(),
			shapedRegion.GetNumberOfPixels(), false );

	// zero-fill
	int zeroFillRO = numberOfRO;
	int zeroFillPE = numberOfPE;

	if ( zeroFill.size() == 2 ) {
		zeroFillRO = zeroFill[0];
		zeroFillPE = zeroFill[1];
	}

	// sort slices
	std::vector< float > pss = fid::Common::GetSlices( fid->GetProcpar() );
	std::vector< int > slices = fid::Common::GetSliceTable( pss );

	// prepare output image
	OutputImageType::Pointer outputImage;
	OutputImageType::SizeType size;
	OutputImageType::IndexType index;
	OutputImageType::RegionType region;

	getOutputImage( zeroFillPE, zeroFillRO, //getOutputImage( numberOfPE, numberOfRO,
			numberOfZ, numberOfEchos, zeroFillPE, zeroFillRO, pss, slices, fid, outputImage, size, index, region );

	// permute: [64, 25, 24, 64] -> [64, 64, 25, 24]
	int permutation[] = { 0, 3, 1, 2 };
	shaped = fid::Common::Permute< ComplexType, 4 >( shaped, permutation );

	// --------------------------------------------------------------------------------------
	// - Fourier-transform																	-
	// --------------------------------------------------------------------------------------

	// for each echo
	for ( int t = 0; t < numberOfEchos; ++t ) {
		// for each slice
		for ( int z = 0; z < numberOfZ; ++z ) {
			// fill slice matrix column-wise
			MatrixType sliceMatrix( numberOfPE, numberOfRO );
			FourierType fft( zeroFillPE, zeroFillRO );

			fft.SetOutputKSpace( outputKSpace ); // output kspace yes/no

			// Matrix index...
			shapedIndex[0] = 0;
			shapedIndex[1] = 0;
			shapedIndex[2] = slices[z];
			shapedIndex[3] = t;

			shapedRegion.SetIndex( shapedIndex );

			// Matrix size...
			shapedSize[0] = numberOfPE; //
			shapedSize[1] = numberOfRO; //
			shapedSize[2] = 1;
			shapedSize[3] = 1;

			shapedRegion.SetSize( shapedSize );

			// copy complex data from shaped into the slice-matrix...
			itk::ImageRegionConstIterator< ShapedDataType > iterator( shaped, shapedRegion );
			ComplexType* sliceMatrixPointer = sliceMatrix.data_block();

			// iterate over slice...
			for ( int i = shapedRegion.GetNumberOfPixels(); i > 0; ++sliceMatrixPointer, ++iterator, --i ) {
				( *sliceMatrixPointer ) = iterator.Value();
			}

			// Fourier-transform
			fft.ifft( sliceMatrix );

			// save to output
			index[0] = 0;
			index[1] = 0;
			index[2] = z;
			index[3] = t;

			size[0] = zeroFillPE;
			size[1] = zeroFillRO;
			size[2] = 1;
			size[3] = 1;

			region.SetIndex( index );
			region.SetSize( size );

			itk::ImageRegionIterator< OutputImageType > it( outputImage, region );

			for ( unsigned int column = 0; column < sliceMatrix.cols(); ++column ) {
				for ( unsigned int row = 0; row < sliceMatrix.rows(); ++row, ++it ) {
					it.Set( sliceMatrix( row, column ) );
				}
			}
		}
	}

	// reorient data
	OutputImageType::Pointer reorientedOutputImage = reorientCoronal( outputImage );

	// calculate magnitude
	itk::Image< float, 4 >::Pointer magnitude = fid::Common::ComplexToMagnitude( reorientedOutputImage );

	// save magnitude
	fid::Common::SaveSeries( magnitude, outputFileName );

	if ( outputRealFileName != "" ) {
		fid::Common::SaveSeries( fid::Common::ComplexToReal( outputImage ), outputRealFileName );
	}

	if ( outputImaginaryFileName != "" ) {
		fid::Common::SaveSeries( fid::Common::ComplexToImaginary( outputImage ), outputImaginaryFileName );
	}

	// t1 mapping
	if ( outputFileNameT1 != "" ) {
		calculateT1( outputFileNameT1, magnitude, fid, maxT1, algorithm, outputA, outputB, outputR2, zeroFillPE, zeroFillRO, pss, slices );
	}

	return EXIT_SUCCESS;
}

