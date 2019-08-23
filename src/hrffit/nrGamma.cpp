//#include "nrGamma.h"
//#include "vnl/vnl_math.h"
//
//namespace nr
//{
//
//Gamma::ScalarType Gamma::gammln( const ScalarType xx )
//{
//  int j;
//  ScalarType x, tmp, y, ser;
//  static const ScalarType cof[ 14 ] = {
//    57.1562356658629235,
//    -59.5979603554754912,
//    14.1360979747417471,
//    -0.491913816097620199,
//    .339946499848118887e-4,
//    .465236289270485756e-4,
//    -.983744753048795646e-4,
//    .158088703224912494e-3,
//    -.210264441724104883e-3,
//    .217439618115212643e-3,
//    -.164318106536763890e-3,
//    .844182239838527433e-4,
//    -.261908384015814087e-4,
//    .368991826595316234e-5 };
//
//  if ( xx <= 0 )
//    {
//    throw( "bad arg in gammln" );
//    }
//
//  y = x = xx;
//  tmp = x + 5.24218750000000000;
//  tmp = ( x + 0.5 ) * vcl_log( tmp ) - tmp;
//  ser = 0.999999999999997092;
//  for ( j = 0; j < 14; j++ )
//    {
//    ser += cof[ j ] / ++y;
//    }
//
//  return tmp + vcl_log( 2.5066282746310005 * ser / x );
//}
//
//Gamma::ScalarType Gamma::gamma( const ScalarType& x, const ScalarType& alpha, const ScalarType& beta )
//{
//  return ( vcl_pow( x / beta, alpha - 1 ) * exp( -1. * ( x / beta ) ) ) / ( beta * exp( gammln( alpha ) ) );
//}
//
//Gamma::ScalarType Gamma::hrf( const ScalarType& x, const ScalarType& f, const ScalarType& alpha1, const ScalarType& beta1, const ScalarType& alpha2, const ScalarType& beta2 )
//{
//  return gamma( x, alpha1, beta1 ) - gamma( x, alpha2, beta2 ) / f;
//}
//
//} // end namespace nr
