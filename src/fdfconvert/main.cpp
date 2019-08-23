#include "fdfConvert.h"
#include "tkdCmdParser.h"

int main( int argc, char ** argv )
{
  std::string path;
  std::string output;
  bool verbose = false;

  tkd::CmdParser p( "fdfconvert", "Convert Varian FDF image folder to 3D/4D image" );

  p.AddArgument( path, "input" )
    ->AddAlias( "i" )
    ->SetDescription( "Input FDF directory" )
    ->SetRequired( true );

  p.AddArgument( output, "output" )
    ->AddAlias( "o" )
    ->SetDescription( "Output image" )
    ->SetRequired( true );

  p.AddArgument( verbose, "verbose" )
    ->AddAlias( "v" )
    ->SetDescription( "Verbose output" );

  if ( !p.Parse( argc, argv ) || output == "" || path == "" )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  fdf::FDFConvert convert;
  convert.SetVerbose( verbose );
  convert.Read( path );
  convert.Write( output );

  return 0;
}
