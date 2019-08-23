#include "tkdCmdParser.h"

#include "graphCommon.h"

/**
 * Degree.
 */
class Degree
{

public:

	typedef double ValueType;

	/**
	 * Run.
	 */
	void run( const std::string& input, const std::string& output,
					float threshold, bool total, bool weighted, bool normalize,
					const std::string& maskImageFileName )
	{
		// get dimensions...
		unsigned int dims = graph::Graph< ValueType >::GetImageDimensions( input );

		if ( dims == 2 )
		{
			process( input, output, threshold, total, weighted, normalize, maskImageFileName );

		} else
		{
			std::cerr << "Number of input (matrix) dimensions should be 2!" << std::endl;
			exit( EXIT_FAILURE );
		}
	}

protected:

	/**
	 * Process.
	 */
	void process( const std::string& input, const std::string& output,
					float threshold, bool total, bool weighted, bool normalize,
					const std::string& maskImageFileName )
	{
		graph::Graph< ValueType >::VectorType ks;

		graph::Graph< ValueType >::Degree( ks, input, threshold, total, weighted, normalize, maskImageFileName, output );

		graph::Graph< ValueType >::WriteVectorToFile( output, ks );

		// print average degree ...
		std::cout << "Degree," << ks.mean() << std::endl;
	}
};

/**
 * Main hierarchy.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( "degree", "Calculate (total) degree for given matrix." );

	std::string inputFileName;
	std::string outputFileName;
	std::string maskImageFileName;

	bool total = true;
	bool weighted = false;
	bool normalize = false;
	float threshold = 0.0;

	p.AddArgument( inputFileName, "input" ) ->AddAlias( "i" ) ->SetInput( "filename" ) ->SetDescription( "Input 2D image: adjacency matrix" ) ->SetRequired(
			true )->SetMinMax( 1, 1 );

	p.AddArgument( outputFileName, "output" ) ->AddAlias( "o" ) ->SetInput( "filename" ) ->SetDescription( "Output degree text file" ) ->SetRequired(
			true ) ->SetMinMax( 1, 1 );

	p.AddArgument( maskImageFileName, "mask" ) ->AddAlias( "m" ) ->SetInput( "filename" ) ->SetDescription( "Mask 3D image: map degree values back to mask" ) ->SetRequired(
			false ) ->SetMinMax( 1, 1 );

	p.AddArgument( threshold, "threshold" ) ->AddAlias( "t" ) ->SetInput( "float" ) ->SetDescription( "Threshold [default: 0.0]" ) ->SetRequired(
			false );

	p.AddArgument( total, "total" ) ->AddAlias( "a" ) ->SetInput( "bool" ) ->SetDescription(
			"Calculate total-degree instead of out-degree [default: true]" ) ->SetRequired( false );

	p.AddArgument( weighted, "weighted" ) ->AddAlias( "w" ) ->SetInput( "bool" ) ->SetDescription(
			"Return strength (weighted degree), [default: false]" ) ->SetRequired( false );

	p.AddArgument( normalize, "normalize" ) ->AddAlias( "n" ) ->SetInput( "normalize" ) ->SetDescription(
			"<bool> standardize (subtract mean and divide by sd), [default: false]" ) ->SetRequired( false );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	Degree degree;

	degree.run( inputFileName, outputFileName, threshold, total, weighted, normalize, maskImageFileName );

	return EXIT_SUCCESS;
}



