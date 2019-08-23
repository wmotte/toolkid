#include "tkdCmdParser.h"
#include "annCommon.h"
#include "graphCommon.h"
#include "annRelabel.h"

//Structure of the command line parameters
struct parameters
{
	std::string inputFileName;
	std::string outputFileName;
	std::string maskFileName;
	double t_begin;
	double t_end;
	double t_interval;
	int c_begin;
	int c_end;
	int c_interval;
	int repeat;
	bool verbose;
};

/**
 * Find optimal clusters and threshold.
 */
class AnnFindOptimalClustering
{
public:

	typedef double PixelType;

	typedef itk::Image< unsigned short, 3 > LabelImageType;

	/**
	 * Run.
	 */
	void Run( const parameters& list )
	{
		// if input is 2D image...
		if ( ann::Ann< PixelType >::GetImageDimensions( list.inputFileName ) == 2 )
		{
			RunThresholdLoop( list );
		}
		else
		{
			std::cerr << "*** ERROR ***: Dimensions of input image is not equal to 2 (matrix)!" << std::endl;
			exit( EXIT_FAILURE );
		}
	}

protected:

	/**
	 * Run over thresholds and calculate modularity for each threshold.
	 */
	void RunThresholdLoop( const parameters& list )
	{
		for( double t = list.t_begin; t <= list.t_end; t += list.t_interval )
		{
			if( list.verbose )
				std::cout << "Running threshold: " << t << std::endl;

			graph::Graph< PixelType >::ImagePointerType input;
			graph::Graph< PixelType >::ImagePointerType output;

			// read input matrix ...
			if( list.verbose )
				std::cout << "Read matrix..." << std::endl;

			graph::Graph< PixelType >::ReadMatrix( input, list.inputFileName );

			// threshold input matrix ...
			if( list.verbose )
				std::cout << "Threshold matrix..." << std::endl;

			graph::Graph< PixelType >::GetMatrix( input, output, 0.0, t, true, true );

			// calculate laplacian ...
			if( list.verbose )
				std::cout << "Calculate laplacian..." << std::endl;

			graph::Graph< PixelType >::MatrixType G;
			graph::Graph< PixelType >::GetMatrix( G, output );
			graph::Graph< PixelType >::MatrixType L;
			graph::Graph< PixelType >::DiagMatrixType D;
			graph::Graph< PixelType >::CalculateLaplacian( L, G, D, false );

			// solve eigen system ...
			if( list.verbose )
				std::cout << "Solve generalized eigensystem..." << std::endl;

			graph::Graph< PixelType >::MatrixType eigenVectors;
			graph::Graph< PixelType >::VectorType eigenValues;
			graph::Graph< PixelType >::GeneralizedEigenCalculation( L, D, eigenVectors, eigenValues );

			if( list.verbose )
				std::cout << "Kmeans++ clustering..." << std::endl;

			// output file name ...
			std::stringstream ss;
			ss << list.outputFileName << "_" << t << ".txt";

			// output vectors ...
			std::vector< double > mods;
			std::vector< double > ents;
			std::vector< int > numbs;

			for( int c = list.c_begin; c <= list.c_end; c += list.c_interval )
			{
				if( list.verbose )
					std::cout << "Clusters total: " << c << " -> ";

				// modularity (average of #repeat runs)...
				PixelType modularity;
				PixelType entropy;
				GetAverageModularityAndEntropy( modularity, entropy, G, eigenVectors, c, list, t );

				if( list.verbose )
					std::cout << "modularity, entropy (" << modularity << ", " << entropy << ")" << std::endl;

				mods.push_back( modularity );
				ents.push_back( entropy );
				numbs.push_back( c );
			}
			if( list.verbose )
				std::cout << ")" << std::endl;

			// write output to file...
			Write( mods, ents, numbs, ss.str() );
		}
	}

	/**
	 * Write modality, entropy and cluster-index to txt file.
	 */
	void Write( const std::vector< double >& mods, const std::vector< double >& ents,
			const std::vector< int > numbs, const std::string& outputFileName )
	{
		std::ofstream out( outputFileName.c_str() );

		if ( out.fail() )
		{
			std::cerr << "*** ERROR ***: Not able to write to: " << outputFileName << "!" << std::endl;
			exit( EXIT_FAILURE );
		}

		if( ( mods.size() != ents.size() ) | ( ents.size() != numbs.size() ) )
		{
			std::cerr << "*** ERROR ***: Did not write to: " << outputFileName << "!" << std::endl;
			std::cerr << "Input vector sizes don't match." << std::endl;
			exit( EXIT_FAILURE );
		}

		for ( unsigned int i = 0; i < numbs.size(); i++ )
			out << numbs[i] << "," << mods[i] << "," << ents[i] << std::endl;

		out.close();
	}

	/**
	 * Repeat clustering and output average modularity and entropy outcome.
	 */
	void GetAverageModularityAndEntropy( PixelType& modularity, PixelType& entropy, const graph::Graph< PixelType >::MatrixType& G,
			const graph::Graph< PixelType >::MatrixType& eigenVectors, unsigned int c,
			const parameters& list, double t )
	{
		std::vector< PixelType > modularities;
		std::vector< vnl_vector< PixelType > > totalClusterings;
		for( int i = 0; i < list.repeat; i++ )
		{
			graph::Graph< PixelType >::VectorType clusters;
			graph::Graph< PixelType >::KMeansPP( eigenVectors, c, clusters );

			// convert vnl matrix to std matrix...
			ann::Ann< PixelType >::MatrixType stdG;
			ann::Ann< PixelType >::VectorType stdClusters;
			ann::Ann< PixelType >::ConvertVNL2STD( G, stdG );
			ann::Ann< PixelType >::ConvertVNL2STD( clusters, stdClusters );

			modularities.push_back( ann::Ann< PixelType >::Modularity( stdClusters, stdG ) );
			totalClusterings.push_back( clusters );
		}

		PixelType sum = std::accumulate( modularities.begin(), modularities.end(), static_cast< PixelType >( 0 ) );
		modularity = sum / static_cast< PixelType >( modularities.size() );

		std::stringstream ss;

		ss << list.outputFileName << "_" << t << "_" << c << ".nii.gz";
		entropy = Entropy( totalClusterings, c, list.maskFileName, ss.str() );
	}

	/**
	 * Discrete entropy.
	 */
	PixelType Entropy( const std::vector< vnl_vector< PixelType > >& matrix,
			unsigned int nLabels, const std::string& maskFileName,
			const std::string& outputFileNameRelabel )
	{
		if( matrix.empty() )
		{
			std::cerr << "*** ERROR ***: could not calculate entropy (matrix empty)!" << std::endl;
			exit( EXIT_FAILURE );
		}

		// [1] Convert each vector into 3D image
		// [2] Combine 3D images into 4D image
		// [3] Relabel 4D image into 4D relabeled image
		// [4] Get Entropy from graph::Grapn< PixelType >::GetTimeSeries( XXX );
		// [5] Entropy( VectorType& entropy, const MatrixType& timeSeries, unsigned int nLabels );
		// [6] Return average entropy.

		// [ 1 ] and [ 2 ] and [ 3 ]...
		graph::Graph< PixelType >::Image4DPointerType image4D;
		graph::Graph< PixelType >::Zero4DImage( image4D, maskFileName, matrix.size() );
		graph::Graph< PixelType >::Image3DPointerType ref3D;

		// construct ref image...
		graph::Graph< PixelType >::FillImage( ref3D, maskFileName, matrix[0] );
		graph::Graph< PixelType >::Insert3DImageIn4DImage( image4D, ref3D, 0 );
		Relabel relabel;
		LabelImageType::Pointer refLabel3D;
		relabel.ConvertImageToLabelImage( refLabel3D, ref3D );

		for( unsigned int i = 1; i < matrix.size(); i++ ) // skip first (=ref)...
		{
			graph::Graph< PixelType >::Image3DPointerType image3D;
			graph::Graph< PixelType >::FillImage( image3D, maskFileName, matrix[i] );

			LabelImageType::Pointer label3D;
			relabel.ConvertImageToLabelImage( label3D, image3D );
			LabelImageType::Pointer relabeled3D;
			relabel.RelabelImage( relabeled3D, label3D, refLabel3D );

			relabel.ConvertLabelImageToImage( image3D, relabeled3D );

			graph::Graph< PixelType >::Insert3DImageIn4DImage( image4D, image3D, i );
		}

		// [ DEBUG ] ... write relabeled file to output for check...
		//graph::Graph< PixelType >::Write4DImage( image4D, outputFileNameRelabel );

		// [ 4 ] ...
		graph::Graph< PixelType >::Image3DPointerType mask3D;
		graph::Graph< PixelType >::Load3DImage( mask3D, maskFileName );
		graph::Graph< PixelType >::MatrixType timeSeries;
		graph::Graph< PixelType >::GetTimeSeries( timeSeries, image4D, mask3D );

		// [ 5 ] ...
		graph::Graph< PixelType >::VectorType entropy;

		graph::Graph< PixelType >::Entropy( entropy, timeSeries, nLabels );

		// [ 6 ] ...
		ann::Ann< PixelType >::VectorType stdEntropy;
		ann::Ann< PixelType >::ConvertVNL2STD( entropy, stdEntropy );
		PixelType sum = std::accumulate( stdEntropy.begin(), stdEntropy.end(), static_cast< PixelType >( 0 ) );
		return sum / static_cast< PixelType >( stdEntropy.size() );
	}

	/**
	 * Return discrete entropy.
	 */
	PixelType Entropy( const ann::Ann< PixelType >::VectorType& a, unsigned int nLabels )
	{
		if ( a.empty() )
		{
			std::cerr << "*** ERROR ***: could not calculate entropy (vector empty)!" << std::endl;
			exit( EXIT_FAILURE );
		}

		PixelType H = 0;
		for( unsigned int i = 1; i <= nLabels; i++ )
		{
			PixelType prob = 0;
			for( unsigned int j = 0; j < a.size(); j++ )
			{
				if( a[j] == i )
					prob += 1.;
			}

			prob /= static_cast< PixelType >( a.size() );

			if( prob != 0 )
				H += prob * std::log( prob );
		}
		return -H;
	}
};

/**
 * Modularity.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( argv[0], "Find optimal clusters and threshold using spectral clustering" );

	parameters list;
	list.t_begin = 30;
	list.t_end = 80;
	list.t_interval = 10;
	list.c_begin = 1;
	list.c_end = 20;
	list.c_interval = 1;
	list.repeat = 10;
	list.outputFileName = "threshold";
	list.verbose = false;

	p.AddArgument( list.inputFileName, "input" )->AddAlias( "i" )->SetDescription( "Input file name (2D matrix)" )->SetRequired( true )->SetMinMax( 1, 1 );
	p.AddArgument( list.maskFileName, "mask" )->AddAlias( "m" )->SetDescription( "Mask file name (3D, needed for entropy)" )->SetRequired( true )->SetMinMax( 1, 1 );
	p.AddArgument( list.outputFileName, "output" )->AddAlias( "o" )->SetDescription( "Output file name root (default: 'threshold')" )->SetRequired( false )->SetMinMax( 1, 1 );

	p.AddArgument( list.t_begin, "threshold-begin" )->AddAlias( "tb" )->SetDescription( "Fixed degree threshold begin (default: 30)" )->SetRequired( false )->SetMinMax( 1, 1 );
	p.AddArgument( list.t_end, "threshold-end" )->AddAlias( "te" )->SetDescription( "Fixed degree threshold end (default: 80)" )->SetRequired( false )->SetMinMax( 1, 1 );
	p.AddArgument( list.t_interval, "threshold-interval" )->AddAlias( "ti" )->SetDescription( "Fixed degree threshold interval (default: 10)" )->SetRequired( false )->SetMinMax( 1, 1 );

	p.AddArgument( list.c_begin, "cluster-begin" )->AddAlias( "cb" )->SetDescription( "Cluster begin (default: 1)" )->SetRequired( false )->SetMinMax( 1, 1 );
	p.AddArgument( list.c_end, "cluster-end" )->AddAlias( "ce" )->SetDescription( "Cluster end (default: 20)" )->SetRequired( false )->SetMinMax( 1, 1 );
	p.AddArgument( list.c_interval, "cluster-interval" )->AddAlias( "ci" )->SetDescription( "Cluster interval (default: 1)" )->SetRequired( false )->SetMinMax( 1, 1 );

	p.AddArgument( list.repeat, "repeat" )->AddAlias( "r" )->SetDescription( "Repeat clustering and report average (default: 10)" )->SetRequired( false )->SetMinMax( 1, 1 );

	p.AddArgument( list.verbose, "verbose" )->AddAlias( "v" )->SetDescription( "Print status to screen (default: false)" )->SetRequired( false )->SetMinMax( 1, 1 );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	AnnFindOptimalClustering annFindOptimalClustering;
	annFindOptimalClustering.Run( list );

	return EXIT_SUCCESS;
}

