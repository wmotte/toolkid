#ifndef __tkdMacro_h__
#define __tkdMacro_h__
#include <sstream>

#define throwMacro( x ) \
{ \
  ::std::stringstream ss; \
  ss << x; \
  throw ::itk::ExceptionObject( __FILE__, __LINE__, ss.str().c_str() ); \
}

#endif /*__tkdMacro_h__*/
