#include "tkdCmdParser.h"

#include "annCommon.h"

/**
 * Output minimal connected components.
 */
class AnnMinConnectedComponents
{
public:

	typedef float PixelType;

	/**
	 * Run.
	 */
	void Run( const std::string& inputFileName,
				const std::string& outputFileName,
				const PixelType threshold,
				const unsigned int minVoxels )
	{
		// read image ...
		ann::Ann< PixelType >::ImagePointerType image;
		ann::Ann< PixelType >::GetImage( image, inputFileName );

		// binarize ...
		ann::Ann< PixelType >::LabelImagePointerType binary;
		ann::Ann< PixelType >::ThresholdToBinaryImage( image, threshold, binary );

		// remove isolated components ...
		ann::Ann< PixelType >::LabelImagePointerType components;
		ann::Ann< PixelType >::RemoveIsolatedComponents( binary, components, minVoxels );

		// and write result to disk ...
		ann::Ann< PixelType >::WriteData( components, outputFileName );
	}
};

/**
 * Minimal connected components filter.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( "ann minimal connected component", "Output minimal connected components" );

	std::string inputFileName;
	std::string outputFileName;

	int minVoxels = 500;
	double threshold = 0.5;

	p.AddArgument( inputFileName, "input" )->AddAlias( "i" )->SetDescription( "Input image file" )->SetRequired( true );
	p.AddArgument( outputFileName, "output" )->AddAlias( "o" )->SetDescription( "Output image file" )->SetRequired( true );
	p.AddArgument( minVoxels, "threshold" )->AddAlias( "t" )->SetDescription( "Binarizing threshold (default: 0.5)" )->SetRequired( false );
	p.AddArgument( minVoxels, "min-voxels" )->AddAlias( "min" )->SetDescription( "Minimal voxels in isolated components (default: 500)" )->SetRequired( false );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	AnnMinConnectedComponents minCC;
	minCC.Run( inputFileName, outputFileName, threshold, minVoxels );

	return EXIT_SUCCESS;
}


