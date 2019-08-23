
#include "tkdCmdParser.h"

#include "graphCommon.h"


/**
 * Threshold.
 */
class Threshold
{
public:

	typedef double PixelType;

public:

	/**
	 * Run.
	 */
	void Run( const std::string& inputFileName, const std::string& outputFileName,
			PixelType threshold, PixelType K, bool fixed, bool weightedOutput )
	{
		graph::Graph< PixelType >::WriterType::Pointer writer =
			graph::Graph< PixelType >::WriterType::New();

		writer->SetFileName( outputFileName.c_str() );

		graph::Graph< PixelType >::ImagePointerType image;
		graph::Graph< PixelType >::ReadMatrix( image, inputFileName );
		graph::Graph< PixelType >::ImagePointerType output;
		graph::Graph< PixelType >::GetMatrix( image, output, threshold, K, fixed, weightedOutput );

		writer->SetInput( output );

		try
		{
			writer->Update();
		} catch ( itk::ExceptionObject& e )
		{
			std::cerr << "Error writing: " << outputFileName << std::endl;
			std::cerr << e.GetDescription() << std::endl;
		}
	}
};

/**
 * Main.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( "theshold", "Convert weighted graph into (binary) graph using threshold or fixed K" );

	std::string inputFileName;
	std::string outputFileName;

	float threshold = 0;
	float k = 10;
	bool fixed = false;
	bool weightedOutput = false;

	p.AddArgument( inputFileName, "input" )->AddAlias( "i" )
				->SetDescription( "Input image: 2D adjacency matrix" ) ->SetRequired( true );
	p.AddArgument( outputFileName, "output" ) ->AddAlias( "o" )
				->SetDescription( "Output image: 2D binary matrix" ) ->SetRequired( true );

	p.AddArgument( threshold, "threshold" )->AddAlias( "t" )->SetDescription(
			"Threshold; only include paths with weight in (threshold, inf] (default: 0.0)" );

	p.AddArgument( k, "k" )->AddAlias( "k" )->SetDescription(
			"Fixed K; output binary matrix with fixed average K (default: 10)" );

	p.AddArgument( fixed, "fixed" )->AddAlias( "f" )->SetDescription(
			"Threshold method: fixed K (default: false)" );

	p.AddArgument( weightedOutput, "weighted-output" )->AddAlias( "w" )->SetDescription(
			"Output weighted matrix (default: false)" );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return -1;
	}

	Threshold thres;
	thres.Run( inputFileName, outputFileName, threshold, k, fixed, weightedOutput );

	return EXIT_SUCCESS;
}

