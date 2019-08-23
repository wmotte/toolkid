#ifndef __nrFitLine_h__
#define __nrFitLine_h__
#include "vnl/vnl_vector.h"

namespace nr
{

namespace fit
{

/** (c) Press et al., Numerical Recipes 3rd ed., p.784 **/
template< class T >
void fitline( const vnl_vector< T >& data, T& a, T& b, int start, int end )
{
  const int numberOfPoints = end - start;

  T sx = 0;
  T sy = 0;
  T st2 = 0;

  for( int i = start; i < end; ++i )
    {
    sx += static_cast< T >( i );
    sy += data( i );
    }

  T sxoss = sx / static_cast< T >( numberOfPoints );

  for( int i = start; i < end; ++i )
    {
    T t = static_cast< T >( i ) - sxoss;
    st2 += t * t;
    a += t * data( i );
    }
  a /= st2;
  b = ( sy - sx * a ) / static_cast< T >( numberOfPoints );
}

template< class T >
void fitline( const vnl_vector< T >& data, T& a, T& b )
{
  fitline< T >( data, a, b, 0, data.size() );
}

} // end namespace fit

} // end namespace nr

#endif /*__nrFitLine_h__*/
