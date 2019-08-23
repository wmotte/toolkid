#ifndef __procparser_h__
#define __procparser_h__
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <vector>
#include <string>
#include <sstream>

class Procparser
{
public:
  typedef std::map< std::string, std::vector< std::string > > MapType;

  bool Parse( const std::string& filename );
  bool Has( const std::string& key ) const;
  int GetSize( const std::string& key ) const;
  int GetSize() const;
  std::vector< std::string > GetKeys() const;
  std::vector< std::string > Get( const std::string& key ) const;
  std::string GetString( const std::string& key, int index = 0 ) const;

  template< class T >
  T GetAs( const std::string& key, int index = 0, T d = T() ) const
  {
    if ( Has( key ) && GetSize( key ) > index )
      {
      std::stringstream ss;
      ss << GetString( key, index );
      ss >> d;
      }
    
    return d;
  }

protected:
  bool Try( std::ifstream& in, const std::string base, const std::string extension );
  bool Open( std::ifstream& in, const std::string base );

  MapType m_Procpar;
};

#endif /*__procparser_h__*/
