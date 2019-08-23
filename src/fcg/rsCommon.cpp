#include "rsCommon.h"

namespace rs
{

template< class PixelType >
PixelType
Common< PixelType >::CorrelationCoefficient( const VectorType& a, const VectorType& b )
{
  // cov( a, b ) / ( std( a ) * std( b ) )

  // cov( a, b ):
  // ( a - mean )' * ( b - mean ) / ( n - 1 )

  VectorType c = a - a.mean();
  VectorType d = b - b.mean();

  PixelType cov = dot_product< PixelType >( c, d ) / static_cast< PixelType >( c.size() - 1 );
  PixelType sumC = 0;
  PixelType sumD = 0;
  for( unsigned int i = 0; i < c.size(); ++i )
    {
    sumC += c( i ) * c( i );
    sumD += d( i ) * d( i );
    }

  PixelType stdC = vcl_sqrt( sumC / static_cast< PixelType >( c.size() - 1 ) );
  PixelType stdD = vcl_sqrt( sumD / static_cast< PixelType >( d.size() - 1 ) );

  return ( cov / ( stdC * stdD ) );
}

template< class PixelType >
PixelType
Common< PixelType >::RtoZ( const PixelType& r )
{
  return vcl_log( ( 1.0 + r ) / ( 1.0 - r ) ) * 0.5;
}

template< class PixelType >
PixelType
Common< PixelType >::ZtoR( const PixelType& z )
{
  return ( 1. - vcl_exp( 2. * z ) ) / ( -1. * vcl_exp( 2. * z ) - 1. );
}

} // end namespace rs

template class rs::Common< float >;
template class rs::Common< double >;
