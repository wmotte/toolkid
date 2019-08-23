#include "tkdCmdParser.h"

#include "annCommon.h"


/**
 * Classify using ANN.
 */
class AnnClassify
{
public:

	typedef float PixelType;

	/**
	 * Run.
	 */
	void Run( const std::string& trainingDataFileName, const std::vector< std::string >& imageFileNames,
			const std::string& maskFileName,
			const std::string& outputFileName,
			unsigned int k,
			double errorBound,
			bool normalize,
			bool addSpatialProperties,
			bool addEuclideanDistance,
			bool addAngles,
			std::vector< int > probabilistic,
			double percentage,
			const std::string& dumpFileName,
			bool prioritySearch )
	{
		// [ 1 ] Get training data into matrix and classes into vector ...
		ann::Ann< PixelType >::MatrixType trainingData; // raw data ...
		ann::Ann< PixelType >::VectorType trainingClasses; // classes ...
		ann::Ann< PixelType >::ReadTrainingData( trainingDataFileName, trainingData, trainingClasses, percentage );

		// [ 2 ] Get query data (images) into matrix ...
		ann::Ann< PixelType >::MatrixType queryData;
		ann::Ann< PixelType >::FillIntensityData( queryData, imageFileNames, maskFileName );

		// x, y, z
		if( addSpatialProperties )
		{
			ann::Ann< PixelType >::FillSpatialData( queryData, maskFileName, maskFileName );
		}

		// euclidean distance ...
		if( addEuclideanDistance )
		{
			ann::Ann< PixelType >::FillEuclideanData( queryData, maskFileName, maskFileName );
		}

		// angles ...
		if ( addAngles )
		{
			ann::Ann< PixelType >::FillAnglesData( queryData, maskFileName, maskFileName );
		}

		if( normalize )
		{
			// [ 3 ] Normalize training data ...
			ann::Ann< PixelType >::Normalize( trainingData );

			// [ 4 ] Normalize query data ...
			ann::Ann< PixelType >::Normalize( queryData );
		}

		// [ 5 ] Query ...
		ann::Ann< PixelType >::MatrixType queryResultMatrix;
		ann::Ann< PixelType >::Query( trainingData, trainingClasses, queryData, queryResultMatrix, k, errorBound, probabilistic, dumpFileName, prioritySearch );

		// [ 6 ] Project to image ...
		ann::Ann< PixelType >::Project2Image( queryResultMatrix, maskFileName, outputFileName, probabilistic );
	}
};

/**
 * ANN Classify.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( "ann classify", "Classify images using ANN based on trained data" );

	std::string trainingData;
	std::vector< std::string > inputFileNames;
	std::string maskFileName;
	std::string outputFileName;
	std::string dumpFileName;
	int k = 1;
	double errorBound = 0.0;
	bool normalize = true;

	bool addSpatialProperties = false;
	bool addEuclideanDistance = false;

	std::vector< int > probabilistic;
	double percentage = 100.0;
	bool prioritySearch = false;
	bool  addAngles = false;

	p.AddArgument( trainingData, "training-data" ) ->AddAlias( "t" ) ->SetDescription( "Training data" ) ->SetRequired(
			true );

	p.AddArgument( inputFileNames, "inputs" ) ->AddAlias( "i" ) ->SetDescription( "Input images" ) ->SetRequired( true );

	p.AddArgument( maskFileName, "mask" ) ->AddAlias( "m" ) ->SetDescription( "Mask image" ) ->SetRequired( true );

	p.AddArgument( outputFileName, "output" ) ->AddAlias( "o" ) ->SetDescription( "Output image" ) ->SetRequired( true );

	p.AddArgument( dumpFileName, "dump-tree" ) ->AddAlias( "d" ) ->SetDescription( "Output file; if specified kd-tree is dumped to file [useful for ann2fig]" ) ->SetRequired( false );

	p.AddArgument( k, "k" ) ->AddAlias( "k" ) ->SetDescription( "K (default: 1)" ) ->SetRequired( false );

	p.AddArgument( errorBound, "error-bound" ) ->AddAlias( "e" ) ->SetDescription( "Error boundaries [default: 0.0]" ) ->SetRequired( false );

	p.AddArgument( probabilistic, "probabilistic" ) ->AddAlias( "p" ) ->SetDescription( "Vector with classes (non-zero) to return as probabilisitic classification (Anbeek, Neuroimage 2005)" )->SetRequired( false )->SetMinMax( 1, 1000 );

	p.AddArgument( percentage, "percentage" ) ->AddAlias( "c" ) ->SetDescription( "Use random percentage of training dataset (to save computational burden); [default: 100]" ) ->SetRequired( false );

	p.AddArgument( normalize, "normalize" ) ->AddAlias( "n" ) ->SetDescription( "Normalize features (default: true)" ) ->SetRequired( false );

	p.AddArgument( prioritySearch, "priority-search" ) ->AddAlias( "r" ) ->SetDescription( "Priority search instead of standard search (default: false)" ) ->SetRequired( false );

	p.AddArgument( addSpatialProperties, "spatial-properties" ) ->AddAlias( "s" ) ->SetDescription( "Include spatial properies as features [default: false]" ) ->SetRequired( false );

	p.AddArgument( addEuclideanDistance, "euclidean-distance" ) ->AddAlias( "euc" ) ->SetDescription( "Add euclidean distance as feature, relative to COG) [default: false]" ) ->SetRequired( false );

	p.AddArgument( addAngles, "angles" ) ->AddAlias( "a" ) ->SetDescription( "Include x,y,z angles, relative to COG) [default: false]" ) ->SetRequired( false );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	AnnClassify annClassify;

	annClassify.Run( trainingData, inputFileNames, maskFileName,
			outputFileName, (unsigned) k, errorBound, normalize, addSpatialProperties, addEuclideanDistance, addAngles, probabilistic, percentage, dumpFileName, prioritySearch );

	return EXIT_SUCCESS;
}

