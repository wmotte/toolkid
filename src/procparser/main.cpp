#include "procparser.h"
#include "tkdCmdParser.h"

int main( int argc, char ** argv )
{
  tkd::CmdParser p( "pp", "Read Varian procpar" );

  std::string base;
  std::string key;
  int index = 0;
  bool isCount = false;

  p.AddArgument( base, "fid" )
    ->AddAlias( "f" )
    ->SetInput( "path" )
    ->SetDescription( "Input FID folder" )
    ->SetRequired( true );

  p.AddArgument( key, "key" )
    ->AddAlias( "k" )
    ->SetInput( "string" )
    ->SetDescription( "Parameter key" );

  p.AddArgument( index, "index" )
    ->AddAlias( "i" )
    ->SetInput( "integer" )
    ->SetDescription( "Index" );

  p.AddArgument( isCount, "size" )
    ->AddAlias( "s" )
    ->SetInput( "boolean" )
    ->SetDescription( "Output number of elements" );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  Procparser procparser;
  procparser.Parse( base );

  if ( key == "" )
    {
    std::vector< std::string > keys = procparser.GetKeys();
    for( std::vector< std::string >::iterator i = keys.begin(); i != keys.end(); ++i )
      {
      std::cout << ( *i ) << std::endl;
      }

    return 0;
    }

  if ( isCount )
    {
    std::cout << procparser.GetSize( key ) << std::endl;
    }
  else
    {
    if ( procparser.Has( key ) )
      {
      if ( index > -1 )
        {
        std::cout << procparser.GetString( key, index ) << std::endl;
        }
      else
        {
        for( int i = 0; i < procparser.GetSize( key ); ++i )
          {
          std::cout << procparser.GetString( key, i ) << std::endl;
          }
        }
      }
    }

  return 0;
}
