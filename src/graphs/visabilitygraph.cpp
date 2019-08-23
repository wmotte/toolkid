#include "tkdCmdParser.h"

#include "graphCommon.h"

/**
 * Date: 16-11-2009
 */
class VisabilityGraph
{

public:

	typedef double ValueType;

	/**
	 * Run.
	 */
	void run( const std::string& input, const std::string& output, bool verbose, bool weighted )
	{
		graph::VisabilityGraph< ValueType >::WriteVisabilityMatrix( input, output, verbose, weighted );
	}
};

/**
 * Main hierarchy.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( "visabilitygraph", "Construct visability graph from given time-series." );

	std::string input;
	std::string output;
	bool verbose = false;
	bool weighted = false;

	p.AddArgument( input, "input" ) ->AddAlias( "i" ) ->SetInput( "filename" ) ->SetDescription( "Time-series input (column format)" ) ->SetRequired(
			true ) -> SetMinMax( 1, 1 );

	p.AddArgument( output, "output" ) ->AddAlias( "o" ) ->SetInput( "filename" ) ->SetDescription( "Output association matrix" ) ->SetRequired(
			true ) -> SetMinMax( 1, 1 );

	p.AddArgument( verbose, "verbose" ) ->AddAlias( "v" ) ->SetInput( "bool" ) ->SetDescription( "Verbose [ default: false ]" ) ->SetRequired( false );

	p.AddArgument( weighted, "weighted" ) ->AddAlias( "w" ) ->SetInput( "bool" ) ->SetDescription( "Use slope as edge weight [ default: false ]" ) ->SetRequired( false );

	if ( !p.Parse( argc, argv ) ) {
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	VisabilityGraph visabilityGraph;

	visabilityGraph.run( input, output, verbose, weighted );

	return EXIT_SUCCESS;
}


