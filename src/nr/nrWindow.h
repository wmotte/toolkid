#ifndef __nrWindow_h__
#define __nrWindow_h__
#include "nrUtility.h"

namespace nr
{

namespace window
{
  template< class T >
  struct SquareWindow
  {
    T operator()( int j, int n )
    {
      return 1.;
    }
  };

  template< class T >
  struct WelchWindow
  {
    T operator()( int j, int n )
    {
      return 1. - utility::sqr< T >( 2. * j / ( n - 1. ) - 1. );
    }
  };
  
} // end namespace window

} // end namespace nr

#endif /*__nrWindow_h__*/
