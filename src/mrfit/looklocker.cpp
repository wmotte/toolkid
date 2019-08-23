#include "looklocker.h"

namespace mapping
{

	/**
	 * Run.
	 */
	void T1MappingLookLocker::Run( const cmdParameters& args )
	{
		TimeSeriesType::Pointer input = GetTimeSeries( args.inputFileName );
		FIDPointerType fid = GetFID( args.fidFileName );

		Process( input, fid, args.outputFileName, args.outputAFileName, args.outputBFileName, args.maxT1 );
	}

	/**
	 * Extra function to use from procfid applications.
	 */
	void T1MappingLookLocker::Process( const TimeSeriesType::Pointer& input, const FIDPointerType& fid, const std::string& outputFileName,
        const std::string& outputAFileName, const std::string& outputBFileName, PixelType maxT1 )
	{
		PixelType ti = fid->Procpar< PixelType > ( "ti" );
		PixelType trimage = fid->Procpar< PixelType > ( "imagetime" );

		SliceTableType sliceTable = fid::Common::GetSliceTable( fid::Common::GetSlices( fid->GetProcpar() ) );

		VolumeType::Pointer output = AllocateVolume( input );
		VolumeType::Pointer outputA = AllocateVolume( input );
		VolumeType::Pointer outputB = AllocateVolume( input );

		unsigned int numberOfSlices = GetNumberOfSlices( input );
		unsigned int numberOfVolumes = GetNumberOfVolumes( input );

		// for each slice ...
		for ( unsigned int z = 0; z < numberOfSlices; z++ )
		{
			VectorType times = GetTRs( sliceTable[z], numberOfSlices, numberOfVolumes, trimage, ti );

			// get single slice time-series...
			TimeSeriesType::Pointer timeSeries = GetSliceTimeSeries( z, input );

			// fit slice...
			mrfit::MRFit::Pointer fit = mrfit::MRFit::New();
			fit->SetInput( timeSeries );
			fit->SetTimes( times );
			fit->FitT1( static_cast< mrfit::MRFit::T1FittingType > ( 8 ), maxT1 );

			// insert fitted map into output volume ...
			AddSliceToVolume( fit->GetMap( mrfit::MRFit::T1_MAP ), output, z );
			AddSliceToVolume( fit->GetMap( mrfit::MRFit::A_MAP ), outputA, z );
			AddSliceToVolume( fit->GetMap( mrfit::MRFit::B_MAP ), outputB, z );
		}

		// write output
		Write( output, outputFileName );
	
    	if( ! outputAFileName.empty() )
            Write( outputA, outputAFileName );
		
        if( ! outputBFileName.empty() )
            Write( outputB, outputBFileName );

	}

	/**
	 * Write output to file.
	 */
	void T1MappingLookLocker::Write( const VolumeType::Pointer& input, const std::string& fileName )
	{
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( fileName );
		writer->SetInput( input );

		try
		{
			writer->Update();
		} catch ( itk::ExceptionObject& e )
		{
			std::cerr << "*** ERROR ***: could not write to: " << fileName << std::endl;
			exit( EXIT_FAILURE );
		}
	}

	/**
	 * Return list of TR's for given slice.
	 */
	VectorType T1MappingLookLocker::GetTRs( unsigned int slice, unsigned int z, unsigned int t, PixelType trimage, PixelType ti )
	{
		/*
		 * MATLAB code procmshotepilc.m
		 * ----------------------------
		 *
		 * trimage = read_procpar( imname, 'imagetime' );
		 * ti = read_procpar( imname, 'ti' );
		 *
		 * for i = 1:ns * nv2;
		 *   titable( i ) = ti + ( i - 1 ) * trimage;
		 * end;
		 *
		 * titable = reshape( titable, ns, nv2 );
		 *
		 * for i = 1:ns;
		 *    titemp = squeeze( titable( pssindex( i ), : ) );
		 *  end;
		 */

		unsigned int total = z * t;

		VectorType vector;
		double queue[total];

		for ( unsigned int i = 0; i < total; i++ )
		{
			double value = ti + i * trimage;
			queue[i] = value;
		}

		typedef boost::multi_array< double, 2 > array;
		typedef boost::const_multi_array_ref< double, 2 > const_array_ref;
		typedef boost::general_storage_order< 2 > storage;

		// reshape ...
		array::size_type ordering[] =
		{ 0, 1 };
		bool ascending[] =
		{ true, true };

		// reshape
		boost::array< array::size_type, 2 > dims =
		{
		{ z, t } };
		const_array_ref matrix( queue, dims, storage( ordering, ascending ) );
		matrix.reshape( dims );

		for ( unsigned int i = 0; i < t; i++ )
			vector.push_back( matrix[slice][i] );

		return vector;
	}

	/**
	 * Allocate 3D input, using the 4D input image information.
	 */
	VolumeType::Pointer T1MappingLookLocker::AllocateVolume( const TimeSeriesType::Pointer& input )
	{
		ExtractVolumeFilterType::Pointer filter = ExtractVolumeFilterType::New();
		filter->SetInput( input );

		TimeSeriesType::RegionType region = input->GetLargestPossibleRegion();
		TimeSeriesType::SizeType size = region.GetSize();
		size[3] = 0;
		region.SetSize( size );

		filter->SetExtractionRegion( region );
		filter->Update();

		VolumeType::Pointer volume = filter->GetOutput();
		volume->FillBuffer( 0 );

		return volume;
	}

	/**
	 * Return 4D input data.
	 */
	TimeSeriesType::Pointer T1MappingLookLocker::GetTimeSeries( const std::string& inputFileName )
	{
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( inputFileName );

		try
		{
			reader->Update();
		} catch ( itk::ExceptionObject& e )
		{
			std::cerr << "*** ERROR ***: could not read: " << inputFileName << std::endl;
			exit( EXIT_FAILURE );
		}

		return reader->GetOutput();
	}

	/**
	 * Return all time-slices from one given slice.
	 */
	TimeSeriesType::Pointer T1MappingLookLocker::GetSliceTimeSeries( const unsigned int slice, TimeSeriesType::Pointer input )
	{
		TimeSeriesType::RegionType inputRegion = input -> GetLargestPossibleRegion();
		TimeSeriesType::SizeType inputSize = inputRegion.GetSize();
		TimeSeriesType::IndexType inputIndex = inputRegion.GetIndex();

		TimeSeriesType::IndexType index;
		TimeSeriesType::SizeType size;
		TimeSeriesType::RegionType region;

		index[0] = inputIndex[0];
		index[1] = inputIndex[1];
		index[2] = slice;
		index[3] = inputIndex[3];

		size[0] = inputSize[0];
		size[1] = inputSize[1];
		size[2] = 1;
		size[3] = inputSize[3];

		region.SetIndex( index );
		region.SetSize( size );

		RegionOfInterestImageFilterType::Pointer filter = RegionOfInterestImageFilterType::New();
		filter->SetRegionOfInterest( region );
		filter->SetInput( input );
		filter->Update();

		return filter->GetOutput();
	}

	/**
	 * Return FID pointer.
	 */
	FIDPointerType T1MappingLookLocker::GetFID( const std::string& inputFileName )
	{
		fid::FIDReader::Pointer reader = fid::FIDReader::New();
		reader->SetFileName( inputFileName );

		try
		{
			reader->Read();
		} catch ( itk::ExceptionObject& e )
		{
			std::cerr << "*** ERROR ***: could not read FID: " << inputFileName << std::endl;
			exit( EXIT_FAILURE );
		}

		return reader->GetFID();
	}

	/**
	 * Add slice to map.
	 */
	void T1MappingLookLocker::AddSliceToVolume( const VolumeType::Pointer& input, VolumeType::Pointer& output, unsigned int slice )
	{
		// slice index...
		PasteImageFilterType::InputImageIndexType index;
		index[0] = 0;
		index[1] = 0;
		index[2] = slice;

		PasteImageFilterType::Pointer filter = PasteImageFilterType::New();
		filter->SetSourceRegion( input->GetLargestPossibleRegion() );
		filter->SetDestinationIndex( index );

		filter->SetSourceImage( input );
		filter->SetDestinationImage( output );
		filter->Update();
		output = filter->GetOutput();
	}

	/**
	 * Return number of Slices.
	 */
	unsigned int T1MappingLookLocker::GetNumberOfSlices( const TimeSeriesType::Pointer& input )
	{
		return input->GetLargestPossibleRegion().GetSize()[2];
	}

	/**
	 * Return number of Volumes.
	 */
	unsigned int T1MappingLookLocker::GetNumberOfVolumes( const TimeSeriesType::Pointer& input )
	{
		return input->GetLargestPossibleRegion().GetSize()[3];
	}

}
