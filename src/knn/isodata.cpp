#include "tkdCmdParser.h"
#include "itkTimeProbe.h"
#include "annCommon.h"
#include "isodata_src/IsoData.h"

/**
 * Fast ISODATA.
 *
 * @author: wim@invivonmr.uu.nl; Image Sciences Institute, UMC Utrecht.
 */
class FastIsoData
{

public:

	typedef float PixelType;

	/**
	 * Run ISODATA (KD-tree).
	 */
	void Run( const std::vector< std::string >& inputFileNames, const std::string& labelFileName, const std::string& outputFileName,
			unsigned int numberOfClusters, unsigned int minNumberOfClusterPoints, unsigned int maxNumberOfIterations, double maxStdev,
			double minRequiredDistance, unsigned int maxPair, bool normalize, bool verbose )
	{
		itk::TimeProbe probe;
		probe.Start();

		// [ 1 ] Read images
		std::vector< ann::Ann< PixelType >::ImagePointerType > images;
		ann::Ann< PixelType >::ReadImages( images, inputFileNames );

		ann::Ann< PixelType >::LabelImagePointerType label;
		ann::Ann< PixelType >::ReadImage( label, labelFileName );

		// [ 2 ] Fill matrix.
		ann::Ann< PixelType >::MatrixType data;
		ann::Ann< PixelType >::FillIntensityData( data, images, label );

		if ( normalize )
		{
			// [ 3 ] Normalize matrix.
			ann::Ann< PixelType >::Normalize( data );
		}

		// [ 4 ] Run clustering.
		std::vector< std::vector< int > > clusters = RunIsoData( data, numberOfClusters, minNumberOfClusterPoints, maxStdev,
				minRequiredDistance, maxPair, maxNumberOfIterations, verbose );

		// [ 5 ] Write to file.
		WriteClusters( clusters, outputFileName, labelFileName );

		probe.Stop();
		if ( verbose )
		{
			std::cout << "Total runtime: " << probe.GetMeanTime() << " s." << std::endl;
		}
	}

protected:

	/**
	 * Write clusters to file.
	 */
	void WriteClusters( const std::vector< std::vector< int > >& clusters, const std::string& outputFileName,
			const std::string& referenceFileName )
	{
		if ( clusters.empty() )
		{
			std::cerr << "*** ERROR ***: no clusters to write to file!" << std::endl;
			exit( EXIT_FAILURE );
		}

		std::vector< PixelType > final;

		// init ... -> to one liner...
		for ( unsigned int i = 0; i < clusters.size(); i++ )
		{
			for ( unsigned int j = 0; j < clusters[i].size(); j++ )
			{
				final.push_back( 0 );
			}
		}

		// add clusters to final list ...
		for ( unsigned int i = 0; i < clusters.size(); i++ )
		{
			for ( unsigned int j = 0; j < clusters[i].size(); j++ )
			{
				final[clusters[i][j]] = i + 1;
			}
		}

		WriteVector2Image( referenceFileName, outputFileName, final );
	}

	/**
	 * Write vector to image, with info from reference.
	 */
	void WriteVector2Image( const std::string& referenceImageFileName, const std::string& outputImageFileName,
			const std::vector< PixelType >& data )
	{
		typedef unsigned char MaskPixelType;
		typedef itk::Image< MaskPixelType, 3 > MaskImageType;
		typedef itk::Image< PixelType, 3 > ImageType;
		typedef itk::ImageMaskSpatialObject< 3 > MaskType;
		typedef itk::ImageSpatialObject< 3, PixelType > SpatialType;
		typedef itk::ImageFileReader< MaskImageType > ReaderType;
		typedef itk::ImageFileWriter< ImageType > WriterType;

		ReaderType::Pointer reader = ReaderType::New();
		reader -> SetFileName( referenceImageFileName.c_str() );
		reader -> Update();

		// mask-image (for dims, etc)
		MaskImageType::Pointer maskImage = reader -> GetOutput();
		MaskType::RegionType region = maskImage -> GetLargestPossibleRegion();
		MaskType::SizeType size = region.GetSize();
		unsigned int voxels = region.GetNumberOfPixels();
		const MaskPixelType* dataMask = maskImage->GetPixelContainer()->GetBufferPointer();

		ImageType::Pointer bufferImage = ImageType::New();
		bufferImage->CopyInformation( maskImage );
		bufferImage->SetRegions( region );
		bufferImage->Allocate();
		bufferImage->FillBuffer( 0 );

		PixelType* buffer = bufferImage -> GetPixelContainer() -> GetBufferPointer();

		// disconnect pipe-lines
		reader = 0;

		// fill data buffer with strengths
		unsigned int count = 0;
		for ( unsigned int i = 0; i < voxels; ++i )
		{
			if ( dataMask[i] != 0 )
			{
				buffer[i] = data[count];
				count++;
			} else
			{
				buffer[i] = 0;
			}
		}

		WriterType::Pointer writer = WriterType::New();

		writer->SetFileName( outputImageFileName.c_str() );
		writer->SetInput( bufferImage );

		try
		{
			writer->Update();
		} catch ( itk::ExceptionObject& e )
		{
			std::cerr << "Error writing: " << outputImageFileName << std::endl;
		}
	}

	/**
	 * Run isodata.
	 */
	std::vector< std::vector< int > > RunIsoData( const ann::Ann< PixelType >::MatrixType& data, int nClusters,
			unsigned int minNumberOfSamplesInCluster, double maxStdev, unsigned int minRequiredDistance, unsigned int maxPair,
			int maxIterations, bool verbose )
	{
		if ( data.empty() )
		{
			std::cerr << "*** ERROR ***: Data matrix is empty!" << std::endl;
			exit( EXIT_FAILURE );
		}

		unsigned int nPts = data[0].size();
		unsigned int dim = data.size();

		ANNpointArray all = kmAllocPts( nPts, dim ); // allocate data points

		for ( unsigned int i = 0; i < nPts; i++ )
		{
			for ( unsigned int d = 0; d < dim; d++ )
			{
				all[i][d] = static_cast< double > ( data[d][i] );
			}
		}

		unsigned int nImages = data.size();
		unsigned int pixels = data[0].size();
		unsigned int nSamples = pixels * nImages;

		// initiate IsoData ...
		IsoData ISODATA = IsoData( pixels, nImages, nClusters, minNumberOfSamplesInCluster );

		// set points ...
		ISODATA.SetPoints( all );

		ISODATA.SampleCenters();

		ISODATA.SamplePoints( nSamples );

		ISODATA.BuildKMfilterCenters();

		for ( int iter = 1; iter <= maxIterations; iter++ )
		{
			if ( iter == maxIterations )
			{
				ISODATA.PreFinalClustering();
			} else if ( iter != 1 && iter != maxIterations )
			{
				ISODATA.SetFilterCenters();
			}
			do
			{
				//STEP 2:
				ISODATA.PutInCluster();

				//STEP 3:
				ISODATA.PostAnalyzeClusters( verbose );

				//STEP 4:
				ISODATA.UpdateCenters();
			} while ( ISODATA.WasDeleted() );

			ISODATA.PutInCluster();

			ISODATA.SetImageCenters();

			//STEP 5:
			ISODATA.CalculateAverageDistances();

			//STEP 6:
			ISODATA.OverallAverageDistances();

			//STEP 7:
			int next_step = 8;
			if ( iter == maxIterations )
			{
				minRequiredDistance = 0;
				next_step = 11;
			}

			else if ( ISODATA.GetNumCenters() <= ( nClusters / 2 ) )
			{
				next_step = 8;
			} else if ( ( iter % 2 == 0 ) || ( ISODATA.GetNumCenters() >= 2 * nClusters ) )
			{
				next_step = 11;
			}

			switch ( next_step )
			{
			case 8:
			{
				//STEP 8:
				ISODATA.CalculateSTDVector();

				//STEP 9:
				ISODATA.CalculateVmax();

				//STEP 10:
				std::vector< int > to_split = ISODATA.ShouldSplit( maxStdev );

				if ( to_split.size() != 0 )
				{
					ISODATA.Split( to_split, verbose );

					if ( iter == maxIterations )
						iter = iter - 1;

					to_split.clear();
					break; //go to step 2
				}
			}

			case 11:
			{
				//STEP 11:
				ISODATA.ComputeCenterDistances();

				//STEP 12:
				vector< PairDistanceNode > to_lump = ISODATA.FindLumpCandidates( minRequiredDistance, maxPair );

				//STEP 13:
				if ( to_lump.size() != 0 )
				{
					ISODATA.Lump( to_lump, verbose );
				}
			}
			}
		}

		if ( verbose )
			ISODATA.GenerateReport( &std::cout );

		// deallocate points ...
		kmDeallocPts( all );

		return ISODATA.GetClusters();
	}
};

/**
 * ISODATA with KD-tree data structure.
 */
int main( int argc, char * argv[] )
{

	// arguments...
	std::vector< std::string > inputs;
	std::string output;
	std::string label;

	int numberOfClusters = 5;
	int minNumberOfClusterPoints = 100;
	int maxNumberOfIter = 100;
	double maxStdev = 0.7;
	double minRequiredDistance = 0.7;
	int maxPair = 2;
	bool normalize = true;
	bool verbose = false;

	tkd::CmdParser parser( argv[0], "Efficient ISODATA (KD-tree data structure)" );

	parser.AddArgument( inputs, "inputs" )->AddAlias( "i" )->SetInput( "<strings>" )->SetDescription( "Input images" )->SetRequired( true )->SetMinMax(
			1, 10000 );

	parser.AddArgument( output, "output" )->AddAlias( "o" )->SetInput( "<string>" )->SetDescription( "Output image" )->SetRequired( true )->SetMinMax(
			1, 1 );

	parser.AddArgument( label, "mask" )->AddAlias( "m" )->SetInput( "<string>" )->SetDescription( "Mask image" )->SetRequired( true )->SetMinMax(
			1, 1 );

	parser.AddArgument( numberOfClusters, "clusters" )->AddAlias( "c" )->SetInput( "<int>" )->SetDescription(
			"Initial number of clusters [default: 5]" )->SetMinMax( 1, 1 );

	parser.AddArgument( minNumberOfClusterPoints, "min-clusters" )->AddAlias( "mc" )->SetInput( "<int>" )->SetDescription(
			"Minimum number of points forming a cluster [default: 100]" )->SetMinMax( 1, 1 );

	parser.AddArgument( maxNumberOfIter, "max-iterations" )->AddAlias( "it" )->SetInput( "<int>" )->SetDescription(
			"Maximum number of iterations [default: 100]" )->SetMinMax( 1, 1 );

	parser.AddArgument( maxStdev, "standard-deviation" )->AddAlias( "sd" )->SetInput( "<double>" )->SetDescription(
			"Maximum standard deviation of points from their cluster center along each axis [default: 0.7]" )->SetMinMax( 1, 1 );

	parser.AddArgument( minRequiredDistance, "min-distance" )->AddAlias( "md" )->SetInput( "<double>" )->SetDescription(
			"Minimum required distance between two cluster centers [default: 0.7]" )->SetMinMax( 1, 1 );

	parser.AddArgument( maxPair, "max-pairs" )->AddAlias( "mp" )->SetInput( "<int>" )->SetDescription(
			"Maximum number of cluster pairs that can be merged per iteration [default: 2]" )->SetMinMax( 1, 1 );

	parser.AddArgument( normalize, "normalize" )->AddAlias( "n" )->SetInput( "<bool>" )->SetDescription(
			"Normalize input images before clustering [default: true]" );

	parser.AddArgument( verbose, "verbose" )->AddAlias( "v" )->SetInput( "<bool>" )->SetDescription(
			"Print status, cluster report and total runtime [default: false]" );

	if ( !parser.Parse( argc, argv ) )
	{
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	FastIsoData isoData = FastIsoData();

	isoData.Run( inputs, label, output, numberOfClusters, minNumberOfClusterPoints, maxNumberOfIter, maxStdev, minRequiredDistance,
			maxPair, normalize, verbose );
}

