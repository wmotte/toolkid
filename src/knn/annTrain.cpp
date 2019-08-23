#include "tkdCmdParser.h"

#include "annCommon.h"

/**
 * Create training dataset ANN.
 */
class AnnTrain
{
public:

	typedef float PixelType;

	/**
	 * Run.
	 */
	void Run( const std::vector< std::string >& inputFileNames,
			const std::string& labelFileName,
			const std::string& outputFileName,
			bool addSpatialCoordinates,
			bool addEuclideanDistance,
			bool addAngles,
			const std::string& maskFileName )
	{
		// [ 1 ] Read images
		std::vector< ann::Ann< PixelType >::ImagePointerType > images;
		ann::Ann< PixelType >::ReadImages( images, inputFileNames );

		ann::Ann< PixelType>::LabelImagePointerType label;
		ann::Ann< PixelType >::ReadImage( label, labelFileName );

		// [ 2 ] Fill feature space.
		ann::Ann< PixelType >::MatrixType results;
		ann::Ann< PixelType >::FillIntensityData( results, images, label );

		// [ 3 ] Add Spatial features if specified.
		if ( addSpatialCoordinates )
			ann::Ann< PixelType >::FillSpatialData( results, labelFileName, maskFileName );

		// [ 4 ] Add euclidean distance.
		if ( addEuclideanDistance )
			ann::Ann< PixelType >::FillEuclideanData( results, labelFileName, maskFileName );

		// [ 5 ] Add angles.
		if ( addAngles )
			ann::Ann< PixelType >::FillAnglesData( results, labelFileName, maskFileName );

		// [ 6 ] Add class labels.
		ann::Ann< PixelType >::FillClassLabels( results, label );

		// [ 7 ] Write data to disk.
		ann::Ann< PixelType >::WriteData( results, outputFileName );
	}
};

/**
 * ANN Train.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( "ann train", "Get training data for ANN Classifier" );

	std::vector< std::string > inputFileNames;
	std::string labelFileName;
	std::string outputFileName;
	std::string maskFileName;

	bool addSpatialCoordinates = false;
	bool addEuclideanDistance = false;
	bool addAngles = false;

	p.AddArgument( inputFileNames, "inputs" ) ->AddAlias( "i" ) ->SetDescription( "Input images" ) ->SetRequired( true );

	p.AddArgument( labelFileName, "label" ) ->AddAlias( "l" ) ->SetDescription( "Label image" ) ->SetRequired( true );

	p.AddArgument( outputFileName, "output" ) ->AddAlias( "o" ) ->SetDescription( "Output text file" ) ->SetRequired( true );

	p.AddArgument( maskFileName, "mask" ) ->AddAlias( "m" ) ->SetDescription( "Mask image" ) ->SetRequired( false );

	p.AddArgument( addSpatialCoordinates, "spatial-coordinates" ) ->AddAlias( "s" ) ->SetDescription( "Add spatial coordinates (relative to COG) as features [default: false]" ) ->SetRequired( false );

	p.AddArgument( addEuclideanDistance, "euclidean-distance" ) ->AddAlias( "euc" ) ->SetDescription( "Add euclidean distance as feature, relative to COG) [default: false]" ) ->SetRequired( false );

	p.AddArgument( addAngles, "angles" ) ->AddAlias( "a" ) ->SetDescription( "Add x,y,z angles, relative to COG) [default: false]" ) ->SetRequired( false );


	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	AnnTrain annTrain;

	annTrain.Run( inputFileNames, labelFileName, outputFileName, addSpatialCoordinates, addEuclideanDistance, addAngles, maskFileName );

	return EXIT_SUCCESS;
}

