#include "itkImage.h"
#include <iostream>
#include <fstream>
#include "vnl/vnl_matrix.h"
#include "vnl/algo/vnl_svd.h"
#include "fidFID.h"
#include "fidFIDReader.h"
#include "fidFourier.h"
#include "itkImageRegionIterator.h"
#include "tkdCmdParser.h"
#include "fidEPITools.h"
#include "fidCommon.h"
#include "../dti/dtifit.h"

int main( int argc, char ** argv ) {
	tkd::CmdParser p( "procmsepi", "Fourier transform multi-shot EPI FID" );

	std::string tablePath;
	std::string inputFileName;
	std::string outputFileName;
	std::string outputRealFileName;
	std::string outputImaginaryFileName;
	std::string outputDTI;
	std::string outputDTIExtension = "nii.gz";
	bool outputKSpace = false;
	std::string referenceFileName;
	std::vector< int > zeroFill;
	bool noOrientImage = false;
	bool mirrorGradientTable = false;

	p.AddArgument( inputFileName, "input" ) ->AddAlias( "i" ) ->SetDescription( "FID input folder" ) ->SetRequired( true );

	p.AddArgument( outputFileName, "output" ) ->AddAlias( "o" ) ->SetDescription( "Output image" ) ->SetRequired( true );

	p.AddArgument( outputRealFileName, "output-real" ) ->AddAlias( "or" ) ->SetDescription( "Output real image" );

	p.AddArgument( outputImaginaryFileName, "output-imag" ) ->AddAlias( "oi" ) ->SetDescription( "Output imaginary image" );

	p.AddArgument( outputKSpace, "k-space" ) ->AddAlias( "k" ) ->SetDescription( "Output k-space (no Fourier transform)" );

	p.AddArgument( referenceFileName, "reference" ) ->AddAlias( "r" ) ->SetDescription( "FID reference folder" ) ->SetRequired( true );

	p.AddArgument( referenceFileName, "reference" ) ->AddAlias( "r" ) ->SetDescription( "FID reference folder" ) ->SetRequired( true );

	p.AddArgument( tablePath, "table" ) ->AddAlias( "t" ) ->SetDescription( "EPI indices table" ) ->SetRequired( false );

	p.AddArgument( outputDTI, "output-dti" ) ->AddAlias( "dti" ) ->SetDescription( "Output filename base for DTI reconstruction" );

	p.AddArgument( outputDTI, "extension-dti" ) ->AddAlias( "e" ) ->SetDescription(
			"Output image extension for DTI reconstruction (default: nii.gz)" );

	p.AddArgument( zeroFill, "zerofill" ) ->AddAlias( "zf" ) ->SetInput( "<int int>" ) ->SetDescription( "Zero-fill matrix" ) ->SetMinMax(
			2, 2 );

	p.AddArgument( noOrientImage, "no-orient" ) ->AddAlias( "no" ) ->SetDescription( "Do not re-orient image after processing" );

	p.AddArgument( mirrorGradientTable, "mirror" ) ->AddAlias( "m" ) ->SetDescription( "Add negative directions to gradient table" ) ->SetRequired(
			false );

	if ( !p.Parse( argc, argv ) ) {
		p.PrintUsage( std::cout );
		return -1;
	}

	typedef fid::FID::ComplexType ComplexType;
	typedef fid::FID::PrecisionType PrecisionType;
	typedef fid::FID::BlockType BlockType;
	typedef fid::FID::DataType DataType;
	typedef vnl_matrix< ComplexType > MatrixType;
	typedef itk::Image< ComplexType, 4 > OutputImageType;

	// reference image

	fid::FIDReader::Pointer referenceReader = fid::FIDReader::New();
	referenceReader->SetFileName( referenceFileName );

	try {
		referenceReader->Read();
	} catch ( itk::ExceptionObject& e ) {
		std::cerr << "Error reading reference FID: " << e.GetDescription() << std::endl;
		return -1;
	}

	fid::FID::Pointer reference = referenceReader->GetFID();
	std::vector< int > referenceSlices = fid::Common::GetSliceTable( fid::Common::GetSlices( reference->GetProcpar() ) );

	// load phase-encoding table
	vnl_matrix< int > table;


	if ( tablePath.empty() ) {
		tablePath = fid::Common::GetAbsolutePath( argv[0], std::string( "tables/" ) + reference->Procpar< std::string > ( "ropat" )
			+ ".txt" );
	}

	if ( !fid::Common::LoadTable( tablePath, table ) ) {
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
	if ( zeroFill.size() == 2 ) {
		zeroFillRO = zeroFill[0];
	}

	// rearrange reference FID

	vnl_matrix< ComplexType > referenceMatrix( numberOfZ, nv * np );

	for ( int i = 0; i < reference->GetNumberOfBlocks(); ++i ) {
		BlockType block = reference->GetBlock( i );

		for ( int j = 0; j < nv; ++j ) {
			for ( int k = 0; k < np; ++k ) {
				referenceMatrix.set_column( j * np + k, block.get_column( table( k, j ) ) );
			}
		}
	}

	// build phasemaps (for every slice: numberOfPE x numberOfRO)

	std::vector< vnl_matrix< PrecisionType > > phaseMaps;

	for ( int z = 0; z < numberOfZ; ++z ) {
		vnl_matrix< PrecisionType > phaseMap( numberOfPE, zeroFillRO );
		phaseMap.fill( itk::NumericTraits< PrecisionType >::Zero );

		// get slice
		vnl_vector< ComplexType > row = referenceMatrix.get_row( z );

		// import as matrix
		vnl_matrix_ref< ComplexType > sliceMatrix1( nv, np, row.data_block() );
		vnl_matrix< ComplexType > sliceMatrix( zeroFillRO, numberOfPE );

		for ( int y = 0; y < numberOfPE; ++y ) {
			fid::Fourier< PrecisionType > fft( zeroFillRO );
			fft.SetOutputKSpace( outputKSpace );
			vnl_vector< ComplexType > row = sliceMatrix1.get_row( y );
			fft.ifft( row );
			sliceMatrix.set_column( y, row );
		}

		// compute weights for fit

		vnl_vector< PrecisionType > weights = fid::EPITools::CalculateWeights( sliceMatrix );

		for ( int y = 0; y < numberOfPE; ++y ) {
			phaseMap.set_row( y, fid::EPITools::BuildPhaseMap( sliceMatrix.get_column( y ), weights ) );
		}

		phaseMaps.push_back( phaseMap );
	}

	// read image FID

	fid::FIDReader::Pointer reader = fid::FIDReader::New();
	reader->SetFileName( inputFileName );

	try {
		reader->Read();
	} catch ( itk::ExceptionObject& e ) {
		std::cerr << "Error reading FID: " << e.GetDescription() << std::endl;
		return -1;
	}

	fid::FID::Pointer fid = reader->GetFID();

	int numberOfT = fid->GetNumberOfBlocks();
	int numberOfShots = fid->GetNumberOfTraces() / numberOfZ;
	numberOfPE *= numberOfShots;

	// zero-fill
	int zeroFillPE = numberOfPE;
	if ( zeroFill.size() == 2 ) {
		zeroFillPE = zeroFill[1];
	}

	// prepare output image

	OutputImageType::Pointer outputImage = OutputImageType::New();
	OutputImageType::RegionType region;

	OutputImageType::SizeType size;
	size[0] = zeroFillPE;
	size[1] = zeroFillRO;
	size[2] = numberOfZ;
	size[3] = numberOfT;

	OutputImageType::IndexType index;
	index[0] = 0;
	index[1] = 0;
	index[2] = 0;
	index[3] = 0;

	region.SetIndex( index );
	region.SetSize( size );
	outputImage->SetRegions( region );
	outputImage->Allocate();

	// sort slices

	std::vector< float > pss = fid::Common::GetSlices( fid->GetProcpar() );
	std::vector< int > slices = fid::Common::GetSliceTable( pss );

	// calculate spacing and field of view

	float lro = fid->Procpar< float > ( "lro" );
	float lpe = fid->Procpar< float > ( "lpe" );

	float sum = 0;

	for ( unsigned int i = 1; i < pss.size(); ++i ) {
		sum += ( pss[slices[i]] - pss[slices[i - 1]] );
	}

	float slcdist = pss.size() == 1 ? .1 : ( sum / static_cast< float > ( pss.size() - 1 ) );

	OutputImageType::SpacingType spacing;
	spacing[0] = lpe / static_cast< float > ( size[0] ) * 10.0;
	spacing[1] = lro / static_cast< float > ( size[1] ) * 10.0;
	spacing[2] = slcdist * 10.0;
	spacing[3] = 1;

	outputImage->SetSpacing( spacing );

	float pss0 = fid->Procpar< float > ( "pss0" ) * 10.0;

	OutputImageType::PointType origin;
	for ( int i = 0; i < 3; ++i ) {
		origin[i] = -0.5 * spacing[i] * static_cast< float > ( size[i] );
	}
	origin[2] += pss0;
	origin[3] = 0;

	outputImage->SetOrigin( origin );

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
	for ( int j = 0; j < nv; ++j ) {
		orderedIndex[1] = j;

		for ( int k = 0; k < np; ++k ) {
			orderedIndex[0] = k;
			dataIndex[0] = table( k, j );

			for ( int t = 0; t < numberOfT; ++t ) {
				orderedIndex[3] = t;
				dataIndex[2] = t;

				for ( int z = 0; z < numberOfTraces; ++z ) {
					orderedIndex[2] = z;
					dataIndex[1] = z;

					orderedData->SetPixel( orderedIndex, data->GetPixel( dataIndex ) );
				}
			}
		}
	}

	// re-order FID: RO x PE x slices x blocks
	typedef itk::Image< ComplexType, 4 > ShapedDataType;

	int dimensions1[] = { numberOfRO, numberOfPE / numberOfShots, numberOfZ, numberOfShots, numberOfT };
	int permutation1[] = { 0, 3, 1, 2, 4 };
	int dimensions2[] = { numberOfRO, numberOfPE, numberOfZ, numberOfT };

	ShapedDataType::Pointer shaped = fid::Common::Reshape< ComplexType, 5, 4 >( fid::Common::Permute< ComplexType, 5 >(
			fid::Common::Reshape< ComplexType, 4, 5 >( orderedData, dimensions1 ), permutation1 ), dimensions2, true );

	ShapedDataType::IndexType shapedIndex;
	shapedIndex[0] = 0;
	shapedIndex[1] = 0;

	// for each block
	for ( int t = 0; t < numberOfT; ++t ) {
		shapedIndex[3] = t;

		// for each slice
		for ( int z = 0; z < numberOfZ; ++z ) {
			shapedIndex[2] = slices[z];

			vnl_matrix_ref< ComplexType > sliceMatrix( numberOfPE, np, &( shaped->GetPixel( shapedIndex ) ) );

			// get slice's phasemap
			const vnl_matrix< PrecisionType >& phaseMap = phaseMaps[referenceSlices[z]];

			// intermediate matrix
			vnl_matrix< ComplexType > intermediate( zeroFillRO, numberOfPE );

			int phaseIndex = 0;
			for ( int y = 0; y < numberOfPE; ++y ) {
				vnl_vector< ComplexType > trace = sliceMatrix.get_row( y );
				fid::Fourier< PrecisionType > fft1( zeroFillRO );
				fft1.SetOutputKSpace( outputKSpace );
				fft1.ifft( trace );

				// phase map
				fid::EPITools::PhaseWithMap( trace, phaseMap.get_row( phaseIndex ) );

				// store
				intermediate.set_column( y, trace );

				if ( ( y + 1 ) % numberOfShots == 0 ) {
					++phaseIndex;
				}
			}

			vnl_matrix< ComplexType > intermediate2( zeroFillRO, zeroFillPE );

			for ( int x = 0; x < zeroFillRO; ++x ) {
				vnl_vector< ComplexType > trace = intermediate.get_row( x );
				fid::Fourier< PrecisionType > fft1( zeroFillPE );
				fft1.SetOutputKSpace( outputKSpace );
				fft1.ifft( trace );

				intermediate2.set_row( x, trace );
			}

			// save to output
			index[2] = z;
			index[3] = t;

			size[2] = 1;
			size[3] = 1;

			region.SetIndex( index );
			region.SetSize( size );

			itk::ImageRegionIterator< OutputImageType > it( outputImage, region );

			for ( int row = 0, column = 0; !it.IsAtEnd(); ++it, ++column ) {
				if ( column >= zeroFillPE ) {
					column = 0;
					++row;
				}

				it.Set( intermediate2( row, column ) );
			}
		}
	}

	// reorient data
	if ( !noOrientImage ) {
		outputImage = fid::Common::ReorientCoronal( outputImage, fid->GetProcpar() );
	}

	// write results

	if ( outputFileName != "" ) {
		fid::Common::SaveSeries( fid::Common::ComplexToMagnitude( outputImage ), outputFileName );
	}

	if ( outputRealFileName != "" ) {
		fid::Common::SaveSeries( fid::Common::ComplexToReal( outputImage ), outputRealFileName );
	}

	if ( outputImaginaryFileName != "" ) {
		fid::Common::SaveSeries( fid::Common::ComplexToImaginary( outputImage ), outputImaginaryFileName );
	}

	if ( outputDTI != "" ) {
		dtifit::DTIFit::Pointer dti = dtifit::DTIFit::New();
		//dti->Run( fid::Common::ComplexToMagnitude( outputImage ), fid->GetProcpar(), outputDTI, outputDTIExtension, mirrorGradientTable );
	}

	return 0;
}
