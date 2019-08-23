
#include "vnl/algo/vnl_sparse_symmetric_eigensystem.h"
#include "tkdCmdParser.h"
#include "vnl/vnl_matrix.h"
#include "itkImage.h"

class Test
{

public:

	/**
	 * Run.
	 */
	void Run( const std::string& inputFileName )
	{
		std::cout << "Test..." << std::endl;

		vnl_sparse_symmetric_eigensystem piet;
		vnl_sparse_matrix< double > M;
		int n = 10;

		piet.CalculateNPairs( M, n );

	}

protected:
};

/**
 * Main.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( "test", "Test" );

	std::string inputFileName;
	std::string outputFileName;

	float lowerThreshold = 0;
	float upperThreshold = 1.0;


	p.AddArgument( inputFileName, "input" ) ->AddAlias( "i" ) ->SetDescription( "Input image: 2D undirected adjacency matrix" ) ->SetRequired( true );

	p.AddArgument( upperThreshold, "threshold" ) ->AddAlias( "thu" ) ->AddAlias( "t" ) ->SetDescription(
			"Threshold; only include paths with weight in (0, threshold] (default: 1.0)" );

	p.AddArgument( lowerThreshold, "threshold-lower" ) ->AddAlias( "thl" ) ->SetDescription(
			"Lower threshold; only include paths with weight > threshold (default: 0)" );

	p.AddArgument( outputFileName, "output" ) ->AddAlias( "o" ) ->SetDescription(
			"Output image: minimum spanning tree with weights" ) ->SetRequired( true );


	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return -1;
	}

	Test test;
	test.Run( inputFileName );

	return EXIT_SUCCESS;
}
