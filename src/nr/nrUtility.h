#ifndef __nrUtility_h__
#define __nrUtility_h__
#include "vnl/vnl_vector.h"

namespace nr
{

namespace utility
{

template< class T >
void swap( T* a, T* b )
{
  T c = *a;
  *a = *b;
  *b = c;
}

template< class T >
void swap( T& a, T& b )
{
  T c = a;
  a = b;
  b = c;
}

template< class T >
T sqr( const T& x )
{
  return ( x * x );
}

template< class T >
T mean( const vnl_vector< T >& data )
{
  T sum = 0;
  const int length = data.size();
  for( int i = 0; i < length; ++i )
    {
    sum += data[ i ];
    }
  return sum / static_cast< T >( length );
}

template< class T >
T variance( const vnl_vector< T >& data, const T& mean )
{
  T ssd = 0;
  const int length = data.size();
  for( int i = 0; i < length; ++i )
    {
    T difference = data[ i ] - mean;
    ssd += sqr< T >( difference );
    }
  return ssd / static_cast< T >( length - 1 );
}

template< class T >
T stdev( const vnl_vector< T >& data, const T& mean )
{
  return vcl_sqrt( variance< T >( data, mean ) );
}

template< class T >
void RealToComplexVector( const vnl_vector< T >& input, vnl_vector< T >& output, int n )
{
  output.set_size( n * 2 );
  output.fill( 0 );

  for( unsigned int i = 0; i < input.size(); ++i )
    {
    output[ i * 2 ] = input[ i ];
    }
}

} // end namespace utility

} // end namespace nr

#endif /*__nrUtility_h__*/
