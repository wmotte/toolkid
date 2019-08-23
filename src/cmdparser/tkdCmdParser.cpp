#include "tkdCmdParser.h"
#include <sstream>
#include <typeinfo>

template< class T >
struct StringParser
{
  void operator()( const std::string& input, T& output )
  {
    std::stringstream ss;
    ss << input;
    ss >> output;
  }
};

template<>
struct StringParser< bool >
{
  void operator()( const std::string& input, bool& output )
  {
    if ( input.length() == 0 )
      {
      output = true;
      return;
      }
    else
      {
      if ( input[ 0 ] == '1' || input[ 0 ] == 'T' || input[ 0 ] == 't' || input[ 0 ] == 'y' || input[ 0 ] == 'Y' )
        {
        output = true;
        }
      else
        {
        output = false;
        }
      }
  }
};

namespace tkd
{

std::string VectorFormat( std::string x, bool vector )
{
	if ( !vector )
		{
		std::stringstream ss;
		ss << "<" << x << ">";
		return ss.str();
		}

	std::stringstream ss;
	ss << "<" << x << "> ...";
	return ss.str();
}

template< class T >
struct ClsToString
{
	std::string operator()( bool vector )
	{
		std::stringstream ss;
		ss << typeid( T ).name();
		return VectorFormat( ss.str(), vector );
	}
};

#define clsToStringMacro( t, text ) \
template<> \
struct ClsToString< t > \
{ \
	std::string operator()( bool vector ) \
	{ \
		return VectorFormat( text, vector ); \
	} \
};

clsToStringMacro( std::string, "string" )
clsToStringMacro( float, "float" )
clsToStringMacro( double, "double" )
clsToStringMacro( int, "integer" )
clsToStringMacro( bool, "boolean" )

CmdParser::ArgumentParserBase::ArgList CmdParser::ArgumentParserBase::GetValues( const MapType& map ) const
{
  MapType::const_iterator ptr = map.find( m_Key );
  if ( ptr == map.end() )
    {
    return ArgList();
    }

  return ptr->second;
}

void CmdParser::ArgumentParserBase::Remap( MapType& map )
{
  for( unsigned int i = 0; i < m_Alias.size(); ++i )
    {
    if ( m_Alias[ i ] == m_Key )
      {
      continue;
      }

    MapType::iterator ptr = map.find( m_Alias[ i ] );
    if ( ptr != map.end() )
      {
      ArgList& list = ptr->second;

      if ( map.find( m_Key ) == map.end() )
        {
        map[ m_Key ] = ArgList();
        }

      for( unsigned int j = 0; j < list.size(); ++j )
        {
        map[ m_Key ].push_back( list[ j ] );
        }

      map.erase( map.find( m_Alias[ i ] ) );
      }
    }
}

CmdParser::ArgumentParserBase* CmdParser::ArgumentParserBase::SetDescription( const std::string& description )
{
  m_Description = description;
  return this;
}

CmdParser::ArgumentParserBase* CmdParser::ArgumentParserBase::SetInput( const std::string& input )
{
  m_Input = input;
  return this;
}

CmdParser::ArgumentParserBase* CmdParser::ArgumentParserBase::AddAlias( const std::string& alias )
{
  m_Alias.push_back( alias );
  return this;
}

CmdParser::ArgumentParserBase* CmdParser::ArgumentParserBase::SetRequired( bool required )
{
  m_Required = required;
  return this;
}

CmdParser::ArgumentParserBase* CmdParser::ArgumentParserBase::SetMinMax( int min, int max )
{
  m_Min = min;
  m_Max = max;
  return this;
}

bool CmdParser::ArgumentParserBase::Validate( MapType& map )
{
  this->Remap( map );

  if ( m_Required )
    {
    if ( map.find( m_Key ) == map.end() )
      {
      return false;
      }
    }

  ArgList list = this->GetValues( map );
  if ( ( m_Min != -1 && static_cast< int >( list.size() ) < m_Min ) ||
       ( m_Max != -1 && static_cast< int >( list.size() ) > m_Max ) )
    {
    if ( list.size() > 0 || m_Required )
      {
      return false;
      }
    }

  return true;
}

CmdParser::ArgumentParserBase::ArgumentParserBase()
: m_Required( false ), m_Min( -1 ), m_Max( -1 )
{
}

CmdParser::ArgumentParserBase::~ArgumentParserBase()
{
}

template< class T >
class ArgumentVectorParser : public CmdParser::ArgumentParserBase
{
public:
  ArgumentVectorParser( std::vector< T >& arg, const std::string& key )
    : m_Arg( arg )
  {
    this->m_Key = key;
  }

  virtual bool Validate( MapType& map )
  {
    if ( !ArgumentParserBase::Validate( map ) )
      {
      return false;
      }

    ArgList list = this->GetValues( map );
    m_Arg.clear();

    for( unsigned int i = 0; i < list.size(); ++i )
      {
      T t = T();
      StringParser< T > sp;
      sp( list[ i ], t );
      m_Arg.push_back( t );
      }

    return true;
  }

private:
  std::vector< T >& m_Arg;
};

template< class T >
class ArgumentSingleParser : public CmdParser::ArgumentParserBase
{
public:
  ArgumentSingleParser( T& arg, const std::string& key )
    : m_Arg( arg )
  {
    this->m_Key = key;
  }

  virtual bool Validate( MapType& map )
  {
    if ( !ArgumentParserBase::Validate( map ) )
      {
      return false;
      }

    ArgList list = this->GetValues( map );
    if ( list.size() > 0 )
      {
      StringParser< T > sp;
      sp( list[ 0 ], m_Arg );
      }
    else if ( map.find( m_Key ) != map.end() )
      {
      StringParser< T > sp;
      sp( "", m_Arg );
      }

    return true;
  }

private:
  T& m_Arg;
};

#define addArgumentImplementMacro( t ) \
  CmdParser::ArgumentParserBase* CmdParser::AddArgument( t& arg, const KeyType& key ) \
  { \
    typedef ArgumentSingleParser< t > ParserType; \
    ParserType* p = new ParserType( arg, key ); \
    p->SetInput( ClsToString< t >()( false ) ); \
    this->m_Parsers.push_back( p ); \
    return p; \
  } \
  CmdParser::ArgumentParserBase* CmdParser::AddArgument( std::vector< t >& arg, const KeyType& key ) \
  { \
    typedef ArgumentVectorParser< t > ParserType; \
    ParserType* p = new ParserType( arg, key ); \
    p->SetInput( ClsToString< t >()( true ) ); \
    this->m_Parsers.push_back( p ); \
    return p; \
  }

addArgumentImplementMacro( std::string )
addArgumentImplementMacro( bool )
addArgumentImplementMacro( int )
addArgumentImplementMacro( float )
addArgumentImplementMacro( double )

CmdParser::CmdParser( const std::string& applicationName, const std::string& applicationDescription )
: m_ApplicationName( applicationName ), m_ApplicationDescription( applicationDescription )
{
}

bool CmdParser::Parse( int argc, char ** argv )
{
  this->Initialize( argc, argv );
  return this->Validate();
}

CmdParser::~CmdParser()
{
  for( std::vector< ArgumentParserBase* >::iterator i = m_Parsers.begin(); i != m_Parsers.end(); ++i )
    {
    ArgumentParserBase* parser = *i;
    delete parser;
    }

  m_Parsers.clear();
}

CmdParser::MapType& CmdParser::GetMap()
{
  return m_Map;
}

const CmdParser::MapType& CmdParser::GetMap() const
{
  return m_Map;
}

CmdParser::ArgList& CmdParser::Get( const KeyType& key )
{
  return m_Map[ key ];
}

const CmdParser::ArgList& CmdParser::Get( const KeyType& key ) const
{
  return m_Map.find( key )->second;
}

bool CmdParser::Has( const KeyType& key ) const
{
  return m_Map.find( key ) != m_Map.end();
}

int CmdParser::Size() const
{
  return static_cast< int >( m_Map.size() );

}

int CmdParser::Size( const KeyType& key ) const
{
  return static_cast< int >( m_Map.find( key )->second.size() );
}

const std::string& CmdParser::GetString( const KeyType& key, int index, const std::string& defaultValue ) const
{
  if ( !this->Has( key ) )
    {
    return defaultValue;
    }

  const ArgList& list = Get( key );
  if ( index >= Size( key ) )
    {
    return defaultValue;
    }

  return list[ index ];
}

void CmdParser::Initialize( int argc, char ** argv )
{
  std::string item;
  std::string key;

  MapType& args = m_Map;
  args.clear();

  for( int i = 1; i < argc; ++i )
    {
    item = argv[ i ];

    if ( item.substr( 0, 1 ) == "-" && item.length() > 1 && ( ( item[ 1 ] >= 'a' && item[ 1 ] <= 'z' ) || ( item[ 1 ] >= 'A' && item[ 1 ] <= 'Z' ) || ( item[ 1 ] == '-' ) ) )
      {
      if ( key.length() > 0 && args.find( key ) == args.end() )
        {
        args[ key ] = std::vector< std::string >();
        }

      if ( item.substr( 0, 2 ) == "--" )
        {
        key = item.substr( 2, item.length() - 2 );
        }
      else
        {
        key = item.substr( 1, item.length() - 1 );
        }
      }
    else
      {
      args[ key ].push_back( item );
      }
    }

  if ( key.length() > 0 && args.find( key ) == args.end() )
    {
    args[ key ] = std::vector< std::string >();
    }
}

void CmdParser::Print( std::ostream& os ) const
{
  for( MapType::const_iterator i = m_Map.begin(); i != m_Map.end(); ++i )
    {
    os << "Key: '" << i->first << "'" << std::endl;

    const ArgList& list = i->second;

    int count = 1;
    for( ArgList::const_iterator j = list.begin(); j != list.end(); ++j )
      {
      os << "\t" << count++ << ":\t'" << ( *j ) << "'" << std::endl;
      }
    }
}

void CmdParser::PrintUsage( std::ostream& os ) const
{
  if ( m_ApplicationName != "" )
    {
    os << m_ApplicationName;

    if ( m_ApplicationDescription != "" )
      {
      os << ": " << std::endl << m_ApplicationDescription << std::endl << std::endl;
      }
    }

  os << "Usage:" << std::endl;
  
  std::vector< std::string > switches;
  for( std::vector< ArgumentParserBase* >::const_iterator i = m_Parsers.begin(); i != m_Parsers.end(); ++i )
    {
    ArgumentParserBase* p = ( *i );
    std::stringstream ss;
    ss << "  " << ( p->m_Required ? "  " : "[ " ) << ( p->m_Key.length() == 1 ? "-" : "--" ) << p->m_Key;
    for( unsigned int j = 0; j < p->m_Alias.size(); ++j )
      {
      ss << ", -" << ( p->m_Alias[ j ].length() > 1 ? "-" : "" ) << p->m_Alias[ j ];
      }
    if ( p->m_Input != "" )
      {
      ss << " " << p->m_Input;
      }
    switches.push_back( ss.str() );
    }
    
  unsigned int maxSize = 0;
  for( unsigned int i = 0; i < switches.size(); ++i )
    {
    maxSize = maxSize > switches[ i ].length() ? maxSize : switches[ i ].length();
    }
  
  int s = 0;
  for( std::vector< ArgumentParserBase* >::const_iterator i = m_Parsers.begin(); i != m_Parsers.end(); ++i, ++s )
    {
    ArgumentParserBase* p = ( *i );

    os << switches[ s ] << ": ";
    for( unsigned int j = switches[ s ].length(); j < maxSize; ++j )
      {
      os << " ";
      }

    os << p->m_Description << ( p->m_Required ? "  " : " ]" ) << std::endl;
    }
  os << std::endl;
}

CmdParser::ArgList CmdParser::GetKeys() const
{
  ArgList keys;
  for( MapType::const_iterator i = m_Map.begin(); i != m_Map.end(); ++i )
    {
    keys.push_back( i->first );
    }

  return keys;
}

bool CmdParser::Validate()
{
  bool valid = true;
  for( std::vector< ArgumentParserBase* >::iterator i = m_Parsers.begin(); i != m_Parsers.end(); ++i )
    {
    ArgumentParserBase* parser = ( *i );
    valid = valid && parser->Validate( this->m_Map );
    }

  return valid;
}


} // end namespace tkd

