#include "annRelabel.h"
#include "tkdCmdParser.h"

/**
 * Minimal connected components filter.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( argv[0], "Relabel label image according reference label image (3D only)" );

	std::string inputFileName;
	std::string refFileName;
	std::string outputFileName;

	p.AddArgument( inputFileName, "input" )->AddAlias( "i" )->SetDescription( "Input file name" )->SetRequired( true )->SetMinMax( 1, 1 );
	p.AddArgument( refFileName, "reference" )->AddAlias( "r" )->SetDescription( "Reference file name" )->SetRequired( true )->SetMinMax( 1, 1 );
	p.AddArgument( outputFileName, "output" )->AddAlias( "o" )->SetDescription( "Output file name" )->SetRequired( true )->SetMinMax( 1, 1 );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	Relabel relabel;
	relabel.Run( inputFileName, refFileName, outputFileName );

	return EXIT_SUCCESS;
}

