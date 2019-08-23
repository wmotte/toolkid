
#include "looklocker.h"

/**
 * Main.
 */
int main(int argc, char ** argv) {

	tkd::CmdParser p("t1map_looklocker", "T1 mapping look locker");

	mapping::cmdParameters args;
	args.maxT1 = 10;

	p.AddArgument( args.inputFileName, "input" )
		->AddAlias( "i" )
		->SetDescription( "Input 4D image" )
		->SetRequired( true );

	p.AddArgument( args.fidFileName, "fid" )
		->AddAlias( "f" )
		->SetDescription( "Path to FID (needed for TR's per slice)" )
		->SetRequired( true );

	p.AddArgument( args.outputFileName, "output" )
		->AddAlias( "o" )
		->SetDescription( "Output T1 map" )
		->SetRequired( true );
	
    p.AddArgument( args.outputAFileName, "output-A" )
		->AddAlias( "oa" )
		->SetDescription( "Output A map" )
		->SetRequired( false );
	
    p.AddArgument( args.outputBFileName, "output-B" )
		->AddAlias( "ob" )
		->SetDescription( "Output B map" )
		->SetRequired( false );

	p.AddArgument( args.maxT1, "maximum-T1" )
		->AddAlias( "max" )
		->SetDescription( "Maximum T1 (default: 10 sec)" );

	if ( ! p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	mapping::T1MappingLookLocker map;
	map.Run( args );

	exit( EXIT_SUCCESS );
}
