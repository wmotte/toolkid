#include "procparser.h"

bool Procparser::Try( std::ifstream& in, const std::string base, const std::string extension )
{
  std::stringstream ss;
  ss << base << extension;

  in.open( ss.str().c_str(), std::ios_base::in );
  if ( in.is_open() )
    {
    return true;
    }

  in.close();
  return false;
}

bool Procparser::Open( std::ifstream& in, const std::string base )
{
  if ( Try( in, base, ".fid/procpar" ) )
    {
    return true;
    }

  if ( Try( in, base, "/procpar" ) )
    {
    return true;
    }

  if ( Try( in, base, "" ) )
    {
    return true;
    }

  return false;
}

bool Procparser::Parse( const std::string& filename )
{
  MapType& procpar = m_Procpar;
  std::ifstream in;
  if ( !Open( in, filename ) )
    {
    std::cerr << "Could not read " << filename << std::endl;
    return false;
    }

  while( !in.eof() )
    {
    std::string id;
    std::string stuff;

    in >> id;

    if ( id == "" )
      {
      break;
      }

    for( int i = 0; i < 10; ++i )
      {
      in >> stuff;
      }

    int numberOfValues = 0;

    for( int j = 0; j < 2; ++j )
      {
      in >> numberOfValues;

      for( int i = 0; i < numberOfValues; ++i )
        {
        std::string value;
        in >> value;

        if ( value[ 0 ] == '"' )
          {
          int delimiters;

          do
            {
            delimiters = 0;

            for( unsigned int k = 0; k < value.length(); ++k )
              {
              if ( value[ k ] == '"' )
                {
                ++delimiters;
                }
              }

            if ( delimiters % 2 != 0 )
              {
              std::string add;
              in >> add;
              value = value + std::string( " " ) + add;
              }
            }
          while( delimiters % 2 != 0 && !in.eof() );
          }

        if ( value.length() > 1 && value.substr( 0, 1 ) == "\"" && value.substr( value.length() - 1, 1 ) == "\"" )
          {
          value = value.substr( 1, value.length() - 2 );
          }

        procpar[ id ].push_back( value );
        }
      }
    }

  in.close();

  return true;
}


bool Procparser::Has( const std::string& key ) const
{
  return m_Procpar.find( key ) != m_Procpar.end();
}


int Procparser::GetSize( const std::string& key ) const
{
  if ( !Has( key ) )
    {
    return 0;
    }

  return m_Procpar.find( key )->second.size();
}

int Procparser::GetSize() const
{
  return m_Procpar.size();
}

std::vector< std::string > Procparser::Get( const std::string& key ) const
{
  return m_Procpar.find( key )->second;
}

std::string Procparser::GetString( const std::string& key, int index ) const
{
  return m_Procpar.find( key )->second[ index ];
}

std::vector< std::string > Procparser::GetKeys() const
{
  std::vector< std::string > keys;
  for( MapType::const_iterator i = m_Procpar.begin(); i != m_Procpar.end(); ++i )
    {
    keys.push_back( i->first );
    }
  return keys;
}
