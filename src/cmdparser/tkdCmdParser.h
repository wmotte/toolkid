#ifndef __tkdCmdParser_h__
#define __tkdCmdParser_h__
#include <map>
#include <vector>
#include <string>
#include <iostream>

namespace tkd
{

#define addArgumentMacro( t ) \
  CmdParser::ArgumentParserBase* AddArgument( t& arg, const KeyType& key ); \
  CmdParser::ArgumentParserBase* AddArgument( std::vector< t >& arg, const KeyType& key );

class CmdParser
{
public:
  typedef std::string KeyType;
  typedef std::vector< KeyType > ArgList;
  typedef std::map< KeyType, ArgList > MapType;

  class ArgumentParserBase
  {
  public:
    typedef std::vector< std::string > ArgList;
    typedef std::map< std::string, ArgList > MapType;

    ArgumentParserBase();
    virtual ~ArgumentParserBase();

    virtual bool Validate( MapType& map );

    ArgumentParserBase* SetDescription( const std::string& description );
    ArgumentParserBase* AddAlias( const std::string& alias );
    ArgumentParserBase* SetRequired( bool required );
    ArgumentParserBase* SetMinMax( int min, int max );
    ArgumentParserBase* SetInput( const std::string& input );

  protected:
    ArgList GetValues( const MapType& map ) const;
    void Remap( MapType& map );

    std::string m_Description;
    std::string m_Key;
    std::string m_Input;
    std::vector< std::string > m_Alias;
    bool m_Required;
    int m_Min;
    int m_Max;

    friend class CmdParser;

  private:
    ArgumentParserBase( const ArgumentParserBase& );
    void operator=( const ArgumentParserBase& );
  };

  CmdParser( const std::string& applicationName = "", const std::string& applicationDescription = "" );
  virtual ~CmdParser();

  bool Parse( int argc, char ** argv );
  void Print( std::ostream& os ) const;
  void PrintUsage( std::ostream& os ) const;

  addArgumentMacro( std::string )
  addArgumentMacro( bool )
  addArgumentMacro( float )
  addArgumentMacro( int )
  addArgumentMacro( double )

  ArgList GetKeys() const;
  MapType& GetMap();
  const MapType& GetMap() const;

  ArgList& Get( const KeyType& key );
  const ArgList& Get( const KeyType& key ) const;

  bool Has( const KeyType& key ) const;
  int Size() const;
  int Size( const KeyType& key ) const;

  const std::string& GetString( const KeyType& key, int index = 0, const std::string& defaultValue = "" ) const;

  template< class T >
  T GetAs( const KeyType& key, int index = 0, T defaultValue = T() )
  {
    std::stringstream ss;
    ss << GetString( key, index, defaultValue );
    ss >> defaultValue;
    return defaultValue;
  }

public:
  CmdParser( const CmdParser& );
  void operator=( const CmdParser& );


protected:
  void Initialize( int argc, char ** argv );
  bool Validate();

  MapType m_Map;
  std::vector< ArgumentParserBase* > m_Parsers;
  std::string m_ApplicationName;
  std::string m_ApplicationDescription;

//private:
//  CmdParser( const CmdParser& );
//  void operator=( const CmdParser& );
};

} // end namespace tkd

#endif /*__tkdCmdParser_h__*/
