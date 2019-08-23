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
#include "itkExtractImageFilter.h"
#include <sstream>

namespace procmsepi
{

	/**
	 * procparseOptions
	 */
	struct procparserOptions
	{
		float pps0;
		float lro;
		float lpe;
		int np;
		int nv;
		int ntraces;
		int nblocks;
		int ns;
		int nsteps;
		std::string seqcon;
		int numberOfRO;
		int numberOfPE;
		int numberOfT;
		int zeroFillRO;
		int zeroFillPE;
	};

	/**
	 * Command line parser options.
	 */
	struct cmdOptions
	{
		std::string tablePath;
		std::string inputFileName;
		std::string outputFileName;
		std::string outputRealFileName;
		std::string outputImaginaryFileName;
		std::string referenceFileName;
		std::vector< int > zeroFill;

		bool outputKSpace;
		bool noOrientImage;
		bool mirrorGradientTable;
		bool forceLinear;
		bool forceCentric;
	};

	enum KTrajectoryType
	{
		KOrderNone = 0, KOrderLinear = 1, KOrderCentric = 2
	};

	/**
	 * Procmsepilc5D class.
	 */
	class Procmsepilc
	{

	public:

		typedef fid::FID::ComplexType ComplexType;
		typedef fid::FID::PrecisionType PrecisionType;
		typedef fid::FID::BlockType BlockType;
		typedef fid::FID::DataType DataType;
		typedef vnl_matrix< ComplexType > MatrixType;
		typedef itk::Image< ComplexType, 4 > OutputImageType;
		typedef itk::Image< ComplexType, 4 > OrderedDataType;
		typedef itk::Image< ComplexType, 5 > ShapedDataType;
		typedef std::vector< vnl_matrix< PrecisionType > > PhaseMapsType;
		typedef vnl_matrix< int > TableType;
		typedef std::vector< int > TrajectoryType;
		typedef std::vector< int > RefSlicesType;

		/**
		 * Run.
		 */
		void Run( const cmdOptions& args )
		{
			// get phase maps and phase-table
			TableType table;
			RefSlicesType referenceSlices;
			PhaseMapsType phaseMaps;
			GetPhaseMaps( args, table, phaseMaps, referenceSlices );

			// get K trajectory
			KTrajectoryType kOrder = GetKOrder( args );

			// read image FID
			fid::FID::Pointer fid = GetFIDReader( args.inputFileName )->GetFID();

			// load procparser options
			procparserOptions options = GetProcparserOptions( args, fid, table );

			// prepare output image
			OutputImageType::Pointer outputImage = GetOutputImage( options, fid );

			std::cout << "OutputImage: " << outputImage->GetLargestPossibleRegion() << std::endl;

			OrderedDataType::Pointer orderedData = GetOrderedData( fid->GetData(), options, table );

			// get trajectory
			TrajectoryType trajectory = GetTrajectory( kOrder, options.nsteps, options.nv );

			// TODO DEBUG


			orderedData->Print( std::cout );
			exit( EXIT_FAILURE );



			// uncompressed processing (e.g. multishpt phMRI EPI)
			Process( trajectory, phaseMaps, referenceSlices, options, orderedData, outputImage, fid, args.outputKSpace );

			// reorient data
			if ( !args.noOrientImage )
				outputImage = fid::Common::ReorientCoronal( outputImage, fid->GetProcpar(), true );

			// write results
			WriteOutput( outputImage, args );
		}

		/**
		 * Uncompressed processing (phMRI data, etc).
		 */
		void Process( const TrajectoryType& trajectory, const PhaseMapsType& phaseMaps,
				const RefSlicesType& referenceSlices,
				const procparserOptions& options,
				const OrderedDataType::Pointer& orderedData,
				OutputImageType::Pointer& outputImage,
				const fid::FID::Pointer& fid,
				bool outputKSpace )

		{
			// uncompressed loop
			if( options.seqcon == "ccnsn" )
			{
				int dimensions[] = { options.numberOfRO, options.numberOfPE / options.nsteps, options.ns, options.nsteps, options.nblocks };
				ShapedDataType::Pointer shaped = fid::Common::Reshape< ComplexType, 4, 5 >( orderedData, dimensions );
				UFFT( shaped, trajectory, phaseMaps, referenceSlices, options, outputImage, fid, outputKSpace );
			}

			/*

			//int dimensions[] = { options.numberOfRO, options.numberOfPE / options.nsteps, options.ns, options.nsteps, options.nblocks };
				ShapedDataType::Pointer shaped = fid::Common::Reshape< ComplexType, 4, 5 >( orderedData, dimensions );

				// FFT
				//UFFT( shaped, trajectory, phaseMaps, referenceSlices, options, outputImage, fid, outputKSpace );
			}
			else if( options.seqcon.compare( "ccncn" ) // compressed loop
			{
				//std::cout << "TODO" << std::endl;
				// CFFT( ... );
			}
			else
			{
				std::cout << "*** ERROR ***: unsupported seqcon: " << options.seqcon << "!" << std::endl;
				exit( EXIT_FAILURE );
			}*/
		}

		/**
		 * Uncompressed FFT.
		 */
		void UFFT( const ShapedDataType::Pointer& shaped, const TrajectoryType& trajectory,
				const PhaseMapsType& phaseMaps, const RefSlicesType& referenceSlices,
				const procparserOptions& options, OutputImageType::Pointer& outputImage,
				const fid::FID::Pointer& fid, bool outputKSpace )
		{
			ShapedDataType::IndexType shapedIndex;
			shapedIndex[0] = 0;
			shapedIndex[1] = 0;

			OutputImageType::RegionType region = outputImage->GetLargestPossibleRegion();
			OutputImageType::IndexType index = region.GetIndex();
			OutputImageType::SizeType size = region.GetSize();

			std::vector< int > slices = fid::Common::GetSliceTable( fid::Common::GetSlices( fid->GetProcpar() ) );

			// for each block
			for ( int t = 0; t < options.nblocks; ++t )
			{
				shapedIndex[4] = t;

				// for each slice
				for ( int z = 0; z < options.ns; ++z )
				{
					shapedIndex[2] = slices[z];

					// get slice's phasemap
					const vnl_matrix< PrecisionType >& phaseMap = phaseMaps[referenceSlices[z]];

					// intermediate matrix
					vnl_matrix< ComplexType > intermediate( options.zeroFillRO, options.numberOfPE );

					for ( int s = 0; s < options.nsteps; ++s )
					{
						shapedIndex[3] = s;

						vnl_matrix_ref< ComplexType > sliceMatrix( options.numberOfPE / options.nsteps, options.np, &( shaped->GetPixel( shapedIndex ) ) );

						for ( int y = 0; y < options.numberOfPE / options.nsteps; ++y )
						{
							vnl_vector< ComplexType > trace = sliceMatrix.get_row( y );
							fid::Fourier< PrecisionType > fft1( options.zeroFillRO );
							fft1.SetOutputKSpace( outputKSpace );
							fft1.ifft( trace );

							// phase map
							fid::EPITools::PhaseWithMap( trace, phaseMap.get_row( y ) );

							// store
							intermediate.set_column( trajectory[s * ( options.numberOfPE / options.nsteps ) + y], trace );
						}
					}

					vnl_matrix< ComplexType > intermediate2( options.zeroFillRO, options.zeroFillPE );

					for ( int x = 0; x < options.zeroFillRO; ++x )
					{
						vnl_vector< ComplexType > trace = intermediate.get_row( x );
						fid::Fourier< PrecisionType > fft1( options.zeroFillPE );
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

					for ( int row = 0, column = 0; !it.IsAtEnd(); ++it, ++column )
					{
						if ( column >= options.zeroFillPE )
						{
							column = 0;
							++row;
						}

						it.Set( intermediate2( row, column ) );
					}
				}
			}
		}


		/**
		 * Return trajectory.
		 */
		TrajectoryType GetTrajectory( const KTrajectoryType& kOrder, int numberOfShots, int nv )
		{
			TrajectoryType trajectory;

			if ( kOrder == KOrderCentric )
			{
				std::vector< int > t1;

				for ( int i = 0; i < numberOfShots; ++i )
					t1.push_back( 1 );

				for ( int i = 0; i < numberOfShots / 2; ++i )
					t1[i] = -1;

				std::vector< int > t2;

				for ( int i = -1 * numberOfShots / 2; i <= numberOfShots / 2 - 1; ++i )
					t2.push_back( i );

				for ( int i = 0; i < numberOfShots; ++i )
				{
					for ( int j = 0; j < nv; ++j )
						trajectory.push_back( t2[i] + 0.5 * static_cast< double > ( t1[i] * j * numberOfShots ) );
				}
			}
			else if ( kOrder == KOrderLinear )
			{
				for ( int i = 0; i < numberOfShots; ++i )
				{
					for ( int j = 0; j < nv; ++j )
						trajectory.push_back( -0.5 * static_cast< double > ( numberOfShots * nv ) + i + numberOfShots * j );
				}
			}

			// sort k-space trajectory
			std::map< int, int > sorter;

			for ( unsigned int i = 0; i < trajectory.size(); ++i )
				sorter[trajectory[i]] = i;

			int itTrajectory = 0;

			for ( std::map< int, int >::iterator i = sorter.begin(); i != sorter.end(); ++i, ++itTrajectory )
				trajectory[i->second] = itTrajectory;

			return trajectory;
		}

		/**
		 * Re-order fid data.
		 */
		OrderedDataType::Pointer GetOrderedData( const DataType::Pointer& data, const procparserOptions& options, const TableType& table )
		{
			OrderedDataType::Pointer orderedData = OrderedDataType::New();
			OrderedDataType::RegionType orderedRegion;
			OrderedDataType::SizeType orderedSize;
			OrderedDataType::IndexType orderedIndex;

			orderedSize[0] = options.np;
			orderedSize[1] = options.nv;
			orderedSize[2] = options.ntraces;
			orderedSize[3] = options.numberOfT;

			orderedRegion.SetSize( orderedSize );
			orderedData->SetRegions( orderedRegion );
			orderedData->Allocate();

			DataType::IndexType dataIndex;

			// re-order points
			for ( int j = 0; j < options.nv; ++j )
			{
				orderedIndex[1] = j;

				for ( int k = 0; k < options.np; ++k )
				{
					orderedIndex[0] = k;
					dataIndex[0] = table( k, j );

					for ( int t = 0; t < options.nblocks; ++t )
					{
						orderedIndex[3] = t;
						dataIndex[2] = t;

						for ( int z = 0; z < options.ntraces; ++z )
						{
							orderedIndex[2] = z;
							dataIndex[1] = z;

							orderedData->SetPixel( orderedIndex, data->GetPixel( dataIndex ) );
						}
					}
				}
			}

			return orderedData;
		}

		/**
		 * Get next file name.
		 */
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

		/**
		 * Prepare output image
		 */
		OutputImageType::Pointer GetOutputImage( const procparserOptions& options, const fid::FID::Pointer& fid )
		{
			// sort slices

			std::vector< float > pss = fid::Common::GetSlices( fid->GetProcpar() );
			std::vector< int > slices = fid::Common::GetSliceTable( pss );

			// calculate spacing and field of view


			OutputImageType::Pointer outputImage = OutputImageType::New();
			OutputImageType::RegionType region;

			OutputImageType::SizeType size;
			size[0] = options.zeroFillRO;
			size[1] = options.zeroFillPE;
			size[2] = options.ns;
			size[3] = options.numberOfT;

			OutputImageType::IndexType index;
			index[0] = 0;
			index[1] = 0;
			index[2] = 0;
			index[3] = 0;

			region.SetIndex( index );
			region.SetSize( size );
			outputImage->SetRegions( region );
			outputImage->Allocate();

			// set spacing

			float sum = 0;

			for ( unsigned int i = 1; i < pss.size(); ++i )
			{
				sum += ( pss[slices[i]] - pss[slices[i - 1]] );
			}

			float slcdist = pss.size() == 1 ? .1 : ( sum / static_cast< float > ( pss.size() - 1 ) );

			OutputImageType::SpacingType spacing;
			spacing[0] = options.lpe / static_cast< float > ( size[0] ) * 10.0;
			spacing[1] = options.lro / static_cast< float > ( size[1] ) * 10.0;
			spacing[2] = slcdist * 10.0;
			spacing[3] = 1;

			outputImage->SetSpacing( spacing );

			OutputImageType::PointType origin;
			for ( int i = 0; i < 3; ++i )
			{
				origin[i] = -0.5 * spacing[i] * static_cast< float > ( size[i] );
			}
			origin[2] += options.pps0;
			origin[3] = 0;

			outputImage->SetOrigin( origin );

			return outputImage;
		}

		/**
		 * Write output.
		 */
		void WriteOutput( const OutputImageType::Pointer& outputImage, const cmdOptions& args )
		{
			if ( args.outputFileName != "" )
				fid::Common::SaveSeries( fid::Common::ComplexToMagnitude( outputImage ), args.outputFileName );

			if ( args.outputRealFileName != "" )
				fid::Common::SaveSeries( fid::Common::ComplexToReal( outputImage ), args.outputRealFileName );

			if ( args.outputImaginaryFileName != "" )
				fid::Common::SaveSeries( fid::Common::ComplexToImaginary( outputImage ), args.outputImaginaryFileName );
		}


		/**
		 * Get FID Reader.
		 */
		fid::FIDReader::Pointer GetFIDReader( const std::string& fidFileName )
		{
			fid::FIDReader::Pointer reader = fid::FIDReader::New();
			reader->SetFileName( fidFileName );

			try
			{
				reader->Read();
			} catch ( itk::ExceptionObject& e )
			{
				std::cerr << "Error reading reference FID: " << e.GetDescription() << std::endl;
				exit( EXIT_FAILURE );
			}
			return reader;
		}

		/**
		 * Get phase-encoding table
		 */
		vnl_matrix< int > GetTable( const cmdOptions& args, const fid::FID::Pointer& fid )
		{
			vnl_matrix< int > table;
			std::string tablePath = fid::Common::GetAbsolutePath( args.tablePath, std::string( "tables/" ) + fid->Procpar< std::string > ( "ropat" ) + ".txt" );

			if ( !fid::Common::LoadTable( tablePath, table ) )
			{
				std::cerr << "Can not open table: " << tablePath << std::endl;
				exit( EXIT_FAILURE );
			}
			return table;
		}

		/**
		 * Get type of K-order
		 */
		KTrajectoryType GetKOrder( const cmdOptions& args )
		{
			fid::FID::Pointer reference = GetFIDReader( args.referenceFileName )->GetFID();

			// k-space trajectory
			KTrajectoryType kOrder = KOrderLinear;
			if ( reference->Procpar< std::string > ( "ky_order" ) == "c" )
			{
				kOrder = KOrderCentric;
			} else if ( reference->Procpar< std::string > ( "ky_order" ) == "l" )
			{
				kOrder = KOrderLinear;
			}

			if ( args.forceLinear )
			{
				kOrder = KOrderLinear;
			}

			if ( args.forceCentric )
			{
				kOrder = KOrderCentric;
			}

			return kOrder;
		}

		/**
		 * Return procparser options like, nt, nsteps, nv2, ns, etc.
		 */
		procparserOptions GetProcparserOptions( const cmdOptions& args, const fid::FID::Pointer& fid, const TableType& table )
		{
			procparserOptions options;

			options.seqcon = fid->Procpar< std::string >( "seqcon" );
			options.pps0 = fid->Procpar< float > ( "pss0" ) * 10.0;
			options.lro = fid->Procpar< float > ( "lro" );
			options.lpe = fid->Procpar< float > ( "lpe" );
			int nv2 = fid->Procpar< int > ( "nv2" );
			options.ns = fid->Procpar< int > ( "ns" );
			options.ntraces = fid->GetNumberOfTraces();
			options.nblocks = fid->GetNumberOfBlocks();

			options.np = table.rows();
			options.nv = table.cols();

			options.numberOfRO = options.np;
			options.numberOfPE = options.nv * options.nsteps;

			// zero-fill
			options.zeroFillPE = options.numberOfPE;
			options.zeroFillRO = pow( 2, ceil( log( options.np ) / log( 2. ) ) );

			if ( args.zeroFill.size() == 2 )
				options.zeroFillPE = args.zeroFill[1];

			if( options.seqcon == "ccncn" )
			{
				options.nsteps = options.ntraces / ( options.ns * nv2 ); // Look-locker
				options.numberOfT = options.nblocks * nv2;
			}
			else if( options.seqcon == "ccnsn" )
			{
				options.nsteps = options.ntraces / options.ns; // fMRI
				options.numberOfT = options.nblocks;
			}
			else
			{
				std::cerr << "*** ERROR ***: unsupported seqcon: " << options.seqcon << "!" << std::endl;
				exit( EXIT_FAILURE );
			}

			return options;
		}

		/**
		 * Return phase maps, table and referenceSlices.
		 */
		void GetPhaseMaps( const cmdOptions& args, TableType& table, PhaseMapsType& phaseMaps, RefSlicesType& referenceSlices )
		{
			// reference image

			fid::FIDReader::Pointer referenceReader = fid::FIDReader::New();
			referenceReader->SetFileName( args.referenceFileName );

			try
			{
				referenceReader->Read();
			} catch ( itk::ExceptionObject& e )
			{
				std::cerr << "Error reading reference FID: " << e.GetDescription() << std::endl;
				exit( EXIT_FAILURE );
			}

			fid::FID::Pointer reference = referenceReader->GetFID();

			// k-space trajectory
			KTrajectoryType kOrder = KOrderLinear;
			if ( reference->Procpar< std::string > ( "ky_order" ) == "c" )
			{
				kOrder = KOrderCentric;
			} else if ( reference->Procpar< std::string > ( "ky_order" ) == "l" )
			{
				kOrder = KOrderLinear;
			}

			if ( args.forceLinear )
			{
				kOrder = KOrderLinear;
			}

			if ( args.forceCentric )
			{
				kOrder = KOrderCentric;
			}

			// load table
			table = GetTable( args, reference );

			int np = table.rows();
			int nv = table.cols();
			int numberOfPE = nv;
			int numberOfZ = reference->Procpar< int > ( "ns" ); // slices

			// zero-fill?
			int zeroFillRO = pow( 2, ceil( log( np ) / log( 2. ) ) );
			;
			if ( args.zeroFill.size() == 2 )
			{
				zeroFillRO = args.zeroFill[0];
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
					fft.SetOutputKSpace( args.outputKSpace );
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

			// reference slices
			referenceSlices = fid::Common::GetSliceTable( fid::Common::GetSlices( reference->GetProcpar() ) );
		}
	}; // end class
}
; // end namespace


/**
 * Main.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( "procmsepilc", "Fourier transform multi-shot EPI linear/centric FID" );

	procmsepi::cmdOptions args;

	args.outputKSpace = false;
	args.noOrientImage = false;
	args.mirrorGradientTable = false;
	args.forceLinear = false;
	args.forceCentric = false;
	args.tablePath = argv[0];

	p.AddArgument( args.inputFileName, "input" ) ->AddAlias( "i" ) ->SetDescription( "FID input folder" ) ->SetRequired( true );

	p.AddArgument( args.outputFileName, "output" ) ->AddAlias( "o" ) ->SetDescription( "Output image" ) ->SetRequired( true );

	p.AddArgument( args.outputRealFileName, "output-real" ) ->AddAlias( "or" ) ->SetDescription( "Output real image" );

	p.AddArgument( args.outputImaginaryFileName, "output-imag" ) ->AddAlias( "oi" ) ->SetDescription( "Output imaginary image" );

	p.AddArgument( args.outputKSpace, "k-space" ) ->AddAlias( "k" ) ->SetDescription( "Output k-space (no Fourier transform)" );

	p.AddArgument( args.referenceFileName, "reference" ) ->AddAlias( "r" ) ->SetDescription( "FID reference folder" ) ->SetRequired( true );

	p.AddArgument( args.referenceFileName, "reference" ) ->AddAlias( "r" ) ->SetDescription( "FID reference folder" ) ->SetRequired( true );

	p.AddArgument( args.zeroFill, "zerofill" ) ->AddAlias( "zf" ) ->SetInput( "<int int>" ) ->SetDescription( "Zero-fill matrix" ) ->SetMinMax(
			2, 2 );

	p.AddArgument( args.noOrientImage, "no-orient" ) ->AddAlias( "no" ) ->SetDescription( "Do not re-orient image after processing" );

	p.AddArgument( args.mirrorGradientTable, "mirror" ) ->AddAlias( "m" ) ->SetDescription( "Add negative directions to gradient table" ) ->SetRequired(
			false );

	p.AddArgument( args.forceLinear, "linear" ) ->SetDescription( "Force linear k-space trajectory processing" );

	p.AddArgument( args.forceCentric, "centric" ) ->SetDescription( "Force centric k-space trajectory processing" );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	procmsepi::Procmsepilc app;
	app.Run( args );
	return( EXIT_SUCCESS );
}

