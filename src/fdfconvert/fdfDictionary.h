#ifndef __fdfDictionary_h__
#define __fdfDictionary_h__
#include <map>
#include <vector>
#include <sstream>

namespace fdf
{

/**
 * \class Dictionary
 * \brief Simple dictionary for storing arrayed values, allowing typed retrieval of data
 */
class Dictionary
  : public std::map< std::string, std::vector< std::string > >
{
public:
  typedef std::map< std::string, std::vector< std::string > > Superclass;
  typedef Superclass::iterator iterator;
  typedef Superclass::const_iterator const_iterator;
  
  /**
   * Get value associated with key as type T.
   * If more than one value is stored with key, position may be set to a value > 0
   * In case position is larger than the number of stored values, or 'key' is not found, the default
   * value 'def' is returned.
   */
  template< class T >
  T Get( const std::string& key, int position = 0, T def = T() )
  {
    iterator i = find( key );
    if ( i == end() )
      {
      return def;
      }
      
    if ( i->second.size() <= position )
      {
      return def;
      }
      
    std::stringstream ss;
    ss << i->second[ position ];
    
    try
      {
      ss >> def;
      }
    catch( ... )
      {
      }
      
    return def;
  }
};

} // end namespace fdf

#endif /*__fdfDictionary_h__*/
