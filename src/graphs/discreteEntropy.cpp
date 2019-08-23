#include "tkdCmdParser.h"

#include "graphCommon.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMinimumMaximumImageFilter.h"
#include "graphCommon.h"

/**
 * Discrete entropy calculation from 4D input image.
 *
 * (Shen, 2010; NeuroImage). "Graph-theory based parcellation of functional subunits in the brain from
 * resting-state fMRI data".
 *
 */
class DiscreteEntropy
{

public:

	typedef double PixelType;

	/**
	 * Run.
	 */
	void Run( const std::string& inputFileName, const std::string& maskFileName,
			const std::string& outputFileName, const std::string& outputMajorityVote )
	{
		if ( graph::Graph< PixelType >::GetImageDimensions( inputFileName ) == 4 )
		{
			if ( !outputFileName.empty() )
				graph::Graph< PixelType >::Entropy( inputFileName, maskFileName, outputFileName );

			if( !outputMajorityVote.empty() )
				graph::Graph< PixelType >::MajorityVote( inputFileName, maskFileName, outputMajorityVote );
		}
		else
		{
			std::cerr << "Input should be 4D!" << std::endl;
			exit( EXIT_FAILURE );
		}
	}
};

/**
 * Main.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( argv[0], "Calculate discrete entropy for 4D input label image" );

	std::string inputFileName;
	std::string maskFileName;
	std::string outputFileName;
	std::string outputMajorityVote;

	p.AddArgument( inputFileName, "input" ) ->AddAlias( "i" ) ->SetInput( "filename" ) ->SetDescription(
			"Input 4D image: aligned 3D segmentations" ) ->SetRequired( true );

	p.AddArgument( maskFileName, "mask" ) ->AddAlias( "m" ) ->SetInput( "filename" ) ->SetDescription( "Mask for voxels of interest" ) ->SetRequired(
			true );

	p.AddArgument( outputFileName, "output" ) ->AddAlias( "o" ) ->SetInput( "filename" ) ->SetDescription( "Output entropy image" ) ->SetRequired(
			false );

	p.AddArgument( outputMajorityVote, "output-majority-vote" ) ->AddAlias( "omv" ) ->SetInput( "filename" ) ->SetDescription( "Output group wise majority vote label image" ) ->SetRequired(
			false );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	DiscreteEntropy entropy;

	entropy.Run( inputFileName, maskFileName, outputFileName, outputMajorityVote );

	return EXIT_SUCCESS;
}

