#include "itkImage.h"
#include <iostream>
#include <fstream>
#include "vnl/vnl_matrix.h"
#include "vnl/algo/vnl_svd.h"

#include "fidFID.h"
#include "fidFIDReader.h"
#include "fidFourier.h"
#include "fidEPITools.h"
#include "fidCommon.h"

#include "itkImageRegionIterator.h"
#include "tkdCmdParser.h"

//#include "../dtifit/dtifit.h"
#include <sstream>
#include "looklocker.h"

std::string GetNextFileName( int counter, std::string baseFileName, std::string extension = "nii.gz", int exp = 1000 )
{
	std::stringstream ss;
	ss << baseFileName << ".";
	while ( counter < exp && exp > 1 )
	{
		ss << "0";
		exp /= 10;
	}
	ss << counter;

	if ( extension != "" )
	{
		ss << "." << extension;
	}

	return ss.str();
}

int main( int argc, char ** argv )
{
	tkd::CmdParser p( "procmsepilcc", "Fourier transform multi-shot EPI linear/centric compressed FID" );

	std::string inputFileName;
	std::string outputFileName;
	std::string outputRealFileName;
	std::string outputImaginaryFileName;
	//std::string outputDTI;
	std::string outputT1;
	std::string outputMapExtension = "nii.gz";
	bool outputKSpace = false;
	std::string referenceFileName;
	std::vector< int > zeroFill;
	bool noOrientImage = false;
	bool mirrorGradientTable = false;
	bool forceLinear = false;
	bool forceCentric = false;

	p.AddArgument( inputFileName, "input" ) ->AddAlias( "i" ) ->SetDescription( "FID input folder" ) ->SetRequired( true );

	p.AddArgument( outputFileName, "output" ) ->AddAlias( "o" ) ->SetDescription( "Output image filename base" ) ->SetRequired( false );

	p.AddArgument( outputRealFileName, "output-real" ) ->AddAlias( "or" ) ->SetDescription( "Output real image filename base" );

	p.AddArgument( outputImaginaryFileName, "output-imag" ) ->AddAlias( "oi" ) ->SetDescription( "Output imaginary image filename base" );

	p.AddArgument( outputKSpace, "k-space" ) ->AddAlias( "k" ) ->SetDescription( "Output k-space (no Fourier transform)" );

	p.AddArgument( referenceFileName, "reference" ) ->AddAlias( "r" ) ->SetDescription( "FID reference folder" ) ->SetRequired( true );

	//p.AddArgument( outputDTI, "output-dti" ) ->AddAlias( "dti" ) ->SetDescription( "Output filename base for DTI reconstruction" );

	p.AddArgument( outputT1, "output-t1" ) ->AddAlias( "t1" ) ->SetDescription( "Output filename base for T1 fitting" );

	p.AddArgument( outputMapExtension, "extension" ) ->AddAlias( "extension-map" ) ->AddAlias( "e" ) ->SetDescription(
			"Output map extension (default: nii.gz)" );

	p.AddArgument( zeroFill, "zerofill" ) ->AddAlias( "zf" ) ->SetInput( "<int int>" ) ->SetDescription( "Zero-fill matrix" ) ->SetMinMax(
			2, 2 );

	p.AddArgument( noOrientImage, "no-orient" ) ->AddAlias( "no" ) ->SetDescription( "Do not re-orient image after processing" );

	p.AddArgument( mirrorGradientTable, "mirror" ) ->AddAlias( "m" ) ->SetDescription( "Add negative directions to gradient table" ) ->SetRequired(
			false );

	p.AddArgument( forceLinear, "linear" ) ->SetDescription( "Force linear k-space trajectory processing" );

	p.AddArgument( forceCentric, "centric" ) ->SetDescription( "Force centric k-space trajectory processing" );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return -1;
	}

	enum KTrajectoryType
	{
		KOrderNone = 0, KOrderLinear = 1, KOrderCentric = 2
	};

	typedef fid::FID::ComplexType ComplexType;
	typedef fid::FID::PrecisionType PrecisionType;
	typedef fid::FID::BlockType BlockType;
	typedef fid::FID::DataType DataType;
	typedef vnl_matrix< ComplexType > MatrixType;
	typedef itk::Image< ComplexType, 4 > OutputImageType;

	// reference image

	fid::FIDReader::Pointer referenceReader = fid::FIDReader::New();
	referenceReader->SetFileName( referenceFileName );

	try
	{
		referenceReader->Read();
	} catch ( itk::ExceptionObject& e )
	{
		std::cerr << "Error reading reference FID: " << e.GetDescription() << std::endl;
		return -1;
	}

	fid::FID::Pointer reference = referenceReader->GetFID();
	std::vector< int > referenceSlices = fid::Common::GetSliceTable( fid::Common::GetSlices( reference->GetProcpar() ) );

	// k-space trajectory
	KTrajectoryType kOrder = KOrderLinear;
	if ( reference->Procpar< std::string > ( "ky_order" ) == "c" )
	{
		kOrder = KOrderCentric;
	} else if ( reference->Procpar< std::string > ( "ky_order" ) == "l" )
	{
		kOrder = KOrderLinear;
	}

	if ( forceLinear )
	{
		kOrder = KOrderLinear;
	}

	if ( forceCentric )
	{
		kOrder = KOrderCentric;
	}

	// load phase-encoding table
	vnl_matrix< int > table;
	std::string tablePath = fid::Common::GetAbsolutePath( argv[0], std::string( "tables/" ) + reference->Procpar< std::string > ( "ropat" )
			+ ".txt" );

	if ( !fid::Common::LoadTable( tablePath, table ) )
	{
		std::cerr << "Can not open table: " << tablePath << std::endl;
		return -1;
	}

	int np = table.rows();
	int nv = table.cols();

	int numberOfRO = np;
	int numberOfPE = nv;
	int numberOfZ = reference->Procpar< int > ( "ns" ); // slices

	// zero-fill?
	int zeroFillRO = pow( 2, ceil( log( np ) / log( 2. ) ) );
	;
	if ( zeroFill.size() == 2 )
	{
		zeroFillRO = zeroFill[0];
	}

	// rearrange reference FID

	vnl_matrix< ComplexType > referenceMatrix( numberOfZ, nv * np );

	for ( int i = 0; i < reference->GetNumberOfBlocks(); ++i )
	{
		BlockType block = reference->GetBlock( i );

		for ( int j = 0; j < nv; ++j )
		{
			for ( int k = 0; k < np; ++k )
			{
				// NOTE: reference #rows = shots x slices, but we only need first shot?
				referenceMatrix.set_column( j * np + k, block.get_column( table( k, j ) ).extract( numberOfZ ) );
			}
		}
	}

	// build phasemaps (for every slice: numberOfPE x numberOfRO)

	std::vector< vnl_matrix< PrecisionType > > phaseMaps;

	for ( int z = 0; z < numberOfZ; ++z )
	{
		vnl_matrix< PrecisionType > phaseMap( numberOfPE, zeroFillRO );
		phaseMap.fill( itk::NumericTraits< PrecisionType >::Zero );

		// get slice
		vnl_vector< ComplexType > row = referenceMatrix.get_row( z );

		// import as matrix
		vnl_matrix_ref< ComplexType > sliceMatrix1( nv, np, row.data_block() );
		vnl_matrix< ComplexType > sliceMatrix( zeroFillRO, numberOfPE );

		for ( int y = 0; y < numberOfPE; ++y )
		{
			fid::Fourier< PrecisionType > fft( zeroFillRO );
			fft.SetOutputKSpace( outputKSpace );
			vnl_vector< ComplexType > row = sliceMatrix1.get_row( y );
			fft.ifft( row );
			sliceMatrix.set_column( y, row );
		}

		// compute weights for fit

		vnl_vector< PrecisionType > weights = fid::EPITools::CalculateWeights( sliceMatrix );

		for ( int y = 0; y < numberOfPE; ++y )
		{
			phaseMap.set_row( y, fid::EPITools::BuildPhaseMap( sliceMatrix.get_column( y ), weights ) );
		}

		phaseMaps.push_back( phaseMap );
	}

	// read image FID

	fid::FIDReader::Pointer reader = fid::FIDReader::New();
	reader->SetFileName( inputFileName );

	try
	{
		reader->Read();
	} catch ( itk::ExceptionObject& e )
	{
		std::cerr << "Error reading FID: " << e.GetDescription() << std::endl;
		return -1;
	}

	fid::FID::Pointer fid = reader->GetFID();

	int numberOfT = fid->GetNumberOfBlocks();
	int numberOfImages = fid->Procpar< int > ( "nv2" );
	int numberOfShots = fid->GetNumberOfTraces() / ( numberOfZ * numberOfImages );
	numberOfPE *= numberOfShots;

	// zero-fill
	int zeroFillPE = numberOfPE;
	if ( zeroFill.size() == 2 )
	{
		zeroFillPE = zeroFill[1];
	}

	// prepare output image

	std::vector< OutputImageType::Pointer > outputImages;
	OutputImageType::RegionType region;
	OutputImageType::IndexType index;
	OutputImageType::SizeType size;

	size[0] = zeroFillPE;
	size[1] = zeroFillRO;
	size[2] = numberOfZ;
	size[3] = numberOfImages;

	index[0] = 0;
	index[1] = 0;
	index[2] = 0;
	index[3] = 0;

	for ( int i = 0; i < numberOfT; ++i )
	{
		OutputImageType::Pointer outputImage = OutputImageType::New();

		region.SetIndex( index );
		region.SetSize( size );
		outputImage->SetRegions( region );
		outputImage->Allocate();

		outputImages.push_back( outputImage );
	}

	// sort slices

	std::vector< float > pss = fid::Common::GetSlices( fid->GetProcpar() );
	std::vector< int > slices = fid::Common::GetSliceTable( pss );

	// calculate spacing and field of view

	float lro = fid->Procpar< float > ( "lro" );
	float lpe = fid->Procpar< float > ( "lpe" );

	float sum = 0;

	for ( unsigned int i = 1; i < pss.size(); ++i )
	{
		sum += ( pss[slices[i]] - pss[slices[i - 1]] );
	}

	float slcdist = pss.size() == 1 ? .1 : ( sum / static_cast< float > ( pss.size() - 1 ) );

	OutputImageType::SpacingType spacing;
	spacing[0] = lpe / static_cast< float > ( size[0] ) * 10.0;
	spacing[1] = lro / static_cast< float > ( size[1] ) * 10.0;
	spacing[2] = slcdist * 10.0;
	spacing[3] = 1;

	for ( int i = 0; i < numberOfT; ++i )
	{
		outputImages[i]->SetSpacing( spacing );
	}

	float pss0 = fid->Procpar< float > ( "pss0" ) * 10.0;

	OutputImageType::PointType origin;
	for ( int i = 0; i < 3; ++i )
	{
		origin[i] = -0.5 * spacing[i] * static_cast< float > ( size[i] );
	}
	origin[2] += pss0;
	origin[3] = 0;

	for ( int i = 0; i < numberOfT; ++i )
	{
		outputImages[i]->SetOrigin( origin );
	}

	// order data points

	typedef itk::Image< ComplexType, 4 > OrderedDataType;
	OrderedDataType::Pointer orderedData = OrderedDataType::New();
	OrderedDataType::RegionType orderedRegion;
	OrderedDataType::SizeType orderedSize;
	OrderedDataType::IndexType orderedIndex;

	int numberOfTraces = fid->GetNumberOfTraces();

	orderedSize[0] = np;
	orderedSize[1] = nv;
	orderedSize[2] = numberOfTraces;
	orderedSize[3] = numberOfT;

	orderedRegion.SetSize( orderedSize );
	orderedData->SetRegions( orderedRegion );
	orderedData->Allocate();

	DataType::Pointer data = fid->GetData();
	DataType::IndexType dataIndex;

	// re-order points
	for ( int j = 0; j < nv; ++j )
	{
		orderedIndex[1] = j;

		for ( int k = 0; k < np; ++k )
		{
			orderedIndex[0] = k;
			dataIndex[0] = table( k, j );

			for ( int t = 0; t < numberOfT; ++t )
			{
				orderedIndex[3] = t;
				dataIndex[2] = t;

				for ( int z = 0; z < numberOfTraces; ++z )
				{
					orderedIndex[2] = z;
					dataIndex[1] = z;

					orderedData->SetPixel( orderedIndex, data->GetPixel( dataIndex ) );
				}
			}
		}
	}

	std::vector< int > trajectory;
	if ( kOrder == KOrderCentric )
	{
		std::vector< int > t1;
		for ( int i = 0; i < numberOfShots; ++i )
		{
			t1.push_back( 1 );
		}

		for ( int i = 0; i < numberOfShots / 2; ++i )
		{
			t1[i] = -1;
		}

		std::vector< int > t2;
		for ( int i = -1 * numberOfShots / 2; i <= numberOfShots / 2 - 1; ++i )
		{
			t2.push_back( i );
		}

		for ( int i = 0; i < numberOfShots; ++i )
		{
			for ( int j = 0; j < nv; ++j )
			{
				trajectory.push_back( t2[i] + 0.5 * static_cast< double > ( t1[i] * j * numberOfShots ) );
			}
		}
	} else if ( kOrder == KOrderLinear )
	{
		for ( int i = 0; i < numberOfShots; ++i )
		{
			for ( int j = 0; j < nv; ++j )
			{
				trajectory.push_back( -0.5 * static_cast< double > ( numberOfShots * nv ) + i + numberOfShots * j );
			}
		}
	}

	// sort k-space trajectory
	std::map< int, int > sorter;
	for ( unsigned int i = 0; i < trajectory.size(); ++i )
	{
		sorter[trajectory[i]] = i;
	}

	int itTrajectory = 0;
	for ( std::map< int, int >::iterator i = sorter.begin(); i != sorter.end(); ++i, ++itTrajectory )
	{
		trajectory[i->second] = itTrajectory;
	}

	// re-order FID: RO x PE x slices x blocks
	typedef itk::Image< ComplexType, 6 > ShapedDataType;
	int dimensions1[] =
	{ numberOfRO, numberOfPE / numberOfShots, numberOfZ, numberOfImages, numberOfShots, numberOfT };
	int permutation[] =
	{ 0, 1, 2, 4, 3, 5 };
	int dimensions2[] =
	{ numberOfRO, numberOfPE / numberOfShots, numberOfZ, numberOfShots, numberOfImages, numberOfT };

	ShapedDataType::Pointer shaped = fid::Common::Reshape< ComplexType, 6, 6 >( fid::Common::Permute< ComplexType, 6 >(
			fid::Common::Reshape< ComplexType, 4, 6 >( orderedData, dimensions1 ), permutation ), dimensions2, true );

	ShapedDataType::IndexType shapedIndex;
	shapedIndex[0] = 0;
	shapedIndex[1] = 0;

	// for each block

	for ( int t = 0; t < numberOfT; ++t )
	{
		shapedIndex[5] = t;

		for ( int m = 0; m < numberOfImages; ++m )
		{
			shapedIndex[4] = m;

			// for each slice
			for ( int z = 0; z < numberOfZ; ++z )
			{
				shapedIndex[2] = slices[z];

				// get slice's phasemap
				const vnl_matrix< PrecisionType >& phaseMap = phaseMaps[referenceSlices[z]];

				// intermediate matrix
				vnl_matrix< ComplexType > intermediate( zeroFillRO, numberOfPE );

				for ( int s = 0; s < numberOfShots; ++s )
				{
					shapedIndex[3] = s;

					vnl_matrix_ref< ComplexType > sliceMatrix( numberOfPE / numberOfShots, np, &( shaped->GetPixel( shapedIndex ) ) );

					for ( int y = 0; y < numberOfPE / numberOfShots; ++y )
					{
						vnl_vector< ComplexType > trace = sliceMatrix.get_row( y );
						fid::Fourier< PrecisionType > fft1( zeroFillRO );
						fft1.SetOutputKSpace( outputKSpace );
						fft1.ifft( trace );

						// phase map
						fid::EPITools::PhaseWithMap( trace, phaseMap.get_row( y ) );

						// store
						intermediate.set_column( trajectory[s * ( numberOfPE / numberOfShots ) + y], trace );
					}
				}

				vnl_matrix< ComplexType > intermediate2( zeroFillRO, zeroFillPE );

				for ( int x = 0; x < zeroFillRO; ++x )
				{
					vnl_vector< ComplexType > trace = intermediate.get_row( x );
					fid::Fourier< PrecisionType > fft1( zeroFillPE );
					fft1.SetOutputKSpace( outputKSpace );
					fft1.ifft( trace );

					intermediate2.set_row( x, trace );
				}

				// save to output
				index[2] = z;
				index[3] = m;

				size[2] = 1;
				size[3] = 1;

				region.SetIndex( index );
				region.SetSize( size );

				itk::ImageRegionIterator< OutputImageType > it( outputImages[t], region );

				for ( int row = 0, column = 0; !it.IsAtEnd(); ++it, ++column )
				{
					if ( column >= zeroFillPE )
					{
						column = 0;
						++row;
					}

					it.Set( intermediate2( row, column ) );
				}
			}
		}
	}

	int exponent = pow( 10, floor( vcl_log( numberOfT ) / vcl_log( 10. ) ) );

	for ( int i = 0; i < numberOfT; ++i )
	{
		OutputImageType::Pointer outputImage = outputImages[i];

		// reorient data
		if ( !noOrientImage )
		{
			outputImage = fid::Common::ReorientCoronal( outputImage, fid->GetProcpar(), true );
		}

		// write results

		if ( outputFileName != "" )
		{
			fid::Common::SaveSeries( fid::Common::ComplexToMagnitude( outputImage ), GetNextFileName( i, outputFileName,
					outputMapExtension, exponent ) );
		}

		if ( outputRealFileName != "" )
		{
			fid::Common::SaveSeries( fid::Common::ComplexToReal( outputImage ), GetNextFileName( i, outputRealFileName, outputMapExtension,
					exponent ) );
		}

		if ( outputImaginaryFileName != "" )
		{
			fid::Common::SaveSeries( fid::Common::ComplexToImaginary( outputImage ), GetNextFileName( i, outputImaginaryFileName,
					outputMapExtension, exponent ) );
		}

		/*
        if ( outputDTI != "" )
		{
			dtifit::DTIFit::Pointer dti = dtifit::DTIFit::New();
			dti->Run( fid::Common::ComplexToMagnitude( outputImage ), fid->GetProcpar(), GetNextFileName( i, outputDTI, "", exponent ),
					outputMapExtension, mirrorGradientTable );
		}*/

		if( outputT1 != "" )
		{
			mapping::TimeSeriesType::Pointer series = fid::Common::ComplexToMagnitude( outputImage );

			std::string outputName = GetNextFileName( i, outputT1, "nii.gz", exponent );

			mapping::T1MappingLookLocker fit;

			fit.Process( series, fid, outputName, "", "", 10 );
		}
	}

	return 0;
}
