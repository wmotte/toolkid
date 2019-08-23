#ifndef __nrCorrelation_h__
#define __nrCorrelation_h__
#include "nrFourier.h"
#include "nrUtility.h"

namespace nr
{

namespace correlation
{

/**
 * Non-normalized zero-lag cross correlation
 * @param a Original vector 1
 * @param b Original vector 2
 * @param meanA Mean of vector 1
 * @param meanB Mean of vector 2
 * @return correlation value
 */
template< class T >
inline T CrossCorrelation( const vnl_vector< T >& a, const vnl_vector< T >& b, const T& meanA, const T& meanB )
{
  const int length = a.size();

  T sum = 0;
  for( int i = 0; i < length; ++i )
    {
    T product = ( a( i ) - meanA ) * ( b( i ) - meanB );
    sum += product;
    }

  return ( sum / static_cast< T >( length - 1 ) );
}

/**
 * Non-normalized zero-lag cross correlation.
 * @param a Original vector 1
 * @param b Original vector 2
 * @return correlation value
 */
template< class T >
inline T CrossCorrelation( const vnl_vector< T >& a, const vnl_vector< T >& b )
{
  const int length = a.size();

  T sum = 0;
  for( int i = 0; i < length; ++i )
    {
    sum += a( i ) * b( i );
    }

  return ( sum / static_cast< T >( length - 1 ) );
}

/**
 * Normalized correlation coefficient
 * @param a Original vector 1
 * @param b Original vector 2
 * @param meanA Mean of vector 1
 * @param meanB Mean of vector 2
 * @param sdA Standard deviation of vector 1
 * @param sdB Standard deviation of vector 2
 * @return Correlation coefficient
 */
template< class T >
inline T NormalizedCrossCorrelation( const vnl_vector< T >& a, const vnl_vector< T >& b, const T& meanA, const T& meanB, const T& sdA, const T& sdB )
{
  return CrossCorrelation< T >( a, b, meanA, meanB ) / ( sdA * sdB );
}

/**
 * Normalized correlation coefficient
 * @param a Original vector 1
 * @param b Original vector 2
 * @param sdA Standard deviation of vector 1
 * @param sdB Standard deviation of vector 2
 * @return Correlation coefficient
 */
template< class T >
inline T NormalizedCrossCorrelation( const vnl_vector< T >& a, const vnl_vector< T >& b, const T& sdA, const T& sdB )
{
  return CrossCorrelation< T >( a, b ) / ( sdA * sdB );
}

/**
 * Normalized correlation coefficient
 * @param a Original vector 1
 * @param b Original vector 2
 * @return Correlation coefficient
 */
template< class T >
inline T NormalizedCrossCorrelation( const vnl_vector< T >& a, const vnl_vector< T >& b )
{
  T meanA = utility::mean< T >( a );
  T meanB = utility::mean< T >( b );
  T sdA = utility::stdev< T >( a, meanA );
  T sdB = utility::stdev< T >( b, meanB );

  return CrossCorrelation< T >( a, b, meanA, meanB ) / ( sdA * sdB );
}

/**
 * Cross correlation for all lags
 * @param a Vector 1 (of length n)
 * @param b Vector 2 (of length n)
 * @param output Output vector of length n
 */
template< class T >
void Correlation( const vnl_vector< T >& a, const vnl_vector< T >& b, vnl_vector< T >& output )
{
  const int n = a.size();
  vnl_vector< T > temp( n );
  for( int i = 0; i < n; ++i )
    {
    output[ i ] = a[ i ];
    temp[ i ] = b[ i ];
    }
  fourier::FFT< T >::realft( output, 1 );
  fourier::FFT< T >::realft( temp, 1 );

  const T no2 = static_cast< T >( n >> 1 );

  for( int i = 2; i < n; i += 2 )
    {
    T tmp = output[ i ];
    output[ i ]     = ( output[ i ] * temp[ i ] + output[ i + 1 ] * temp[ i + 1 ] ) / no2;
    output[ i + 1 ] = ( output[ i + 1 ] * temp[ i ] - tmp * temp[ i + 1 ] ) / no2;
    }

  output[ 0 ] = output[ 0 ] * temp[ 0 ] / no2;
  output[ 1 ] = output[ 1 ] * temp[ 1 ] / no2;

  fourier::FFT< T >::realft( output, -1 );
}

/**
 * Functional distance (directly related to pearson correlation coefficient (distance(x,y = 2 - 2 * corr(x,y)).
 *
 * The weight is defined using a Gaussian kernel (scale parameter sigma).
 *
 * Reference: IEEE Trans Pat Ana Machine Intel, Vol 22, No 8, August 2000 (Shi and Malik).
 */
template< class T >
inline T FunctionalDistance( const vnl_vector< T >& a, const vnl_vector< T >& b )
{
	const unsigned int n = a.size();
	vnl_vector< T > functionalDistance = vnl_vector< T >( n );
	for( unsigned i = 0; i < n; i++ )
		functionalDistance[ i ] = std::pow( a[ i ] - b[ i ], static_cast< T >( 2 ) );

	T absDistance = 0;
	for( unsigned i = 0; i < n; i++ )
		absDistance += functionalDistance[ i ];

	return std::sqrt( absDistance );
}

/**
 * Correlation coefficient and lag (in samples) where cross correlation is maximal
 * @param a Vector 1
 * @param b Vector 2
 * @param r Output correlation coefficient
 * @param lag Output lag (in samples)
 */
template< class T >
void TimeLaggedNormalizedCrossCorrelation( const vnl_vector< T >& a, const vnl_vector< T >& b, T& r, int& lag )
{
  const int length = a.size();

  T meanA = utility::mean< T >( a );
  T meanB = utility::mean< T >( b );

  T sdA = utility::stdev< T >( a, meanA );
  T sdB = utility::stdev< T >( b, meanB );

  const int length2 = pow( 2, vcl_ceil( vcl_log( length ) / vcl_log( 2 ) ) );

  vnl_vector< T > tmpA( length2 );
  vnl_vector< T > tmpB( length2 );
  vnl_vector< T > tmpO( length2 );

  for( int i = 0; i < length; ++i )
    {
    tmpA[ i ] = a[ i ] - meanA;
    tmpB[ i ] = b[ i ] - meanB;
    }

  for( int i = length; i < length2; ++i )
    {
    tmpA[ i ] = 0;
    tmpB[ i ] = 0;
    }

  Correlation( tmpA, tmpB, tmpO );

  lag = 0;
  r = 0;

  for( int i = 0; i < length2; ++i )
    {
    if ( i == 0 || tmpO[ i ] > r )
      {
      r = tmpO[ i ];
      lag = i;
      }
    }

  r = ( r / static_cast< T >( length - 1 ) ) / ( sdA * sdB );
}

/**
 * Correlation coefficient and lag (in samples) where cross correlation is maximal
 * @param a Vector 1
 * @param b Vector 2
 * @param r Output correlation coefficient
 * @param lag Output lag (in samples)
 */
template< class T >
void TimeLaggedNormalizedCrossCorrelation( const vnl_vector< T >& a, const vnl_vector< T >& b, const T& meanA, const T& meanB, const T& sdA, const T& sdB, T& r, int& lag )
{
  const int length = a.size();

  const int length2 = pow( 2, vcl_ceil( vcl_log( length ) / vcl_log( 2 ) ) );

  vnl_vector< T > tmpA( length2 );
  vnl_vector< T > tmpB( length2 );
  vnl_vector< T > tmpO( length2 );

  for( int i = 0; i < length; ++i )
    {
    tmpA[ i ] = a[ i ] - meanA;
    tmpB[ i ] = b[ i ] - meanB;
    }

  for( int i = length; i < length2; ++i )
    {
    tmpA[ i ] = 0;
    tmpB[ i ] = 0;
    }

  Correlation( tmpA, tmpB, tmpO );

  lag = 0;
  r = 0;

  for( int i = 0; i < length2; ++i )
    {
    if ( i == 0 || tmpO[ i ] > r )
      {
      r = tmpO[ i ];
      lag = i;
      }
    }

  r = ( r / static_cast< T >( length - 1 ) ) / ( sdA * sdB );
}

/**
 * Hilbert transform, returns analytic signal of an input vector. Complex vector input/output stores
 * real and imaginary parts interleaved.
 * @param input Complex vector of length n*2 (with n a power of 2)
 * @param n Length of the vector
 * @param output Hilbert transform
 */
template< class T >
void Hilbert( const vnl_vector< T >& input, const int& n, vnl_vector< T >& output )
{
  const int n2 = 2 * n;
  output = input;
  nr::fourier::FFT< T >::four1( output.data_block(), n, 1 );

  // k=1 ... n/2
  for( int k = 2; k < n; ++k )
    {
    output( k ) *= 2.;
    }

  // k = n/2+1 ... n
  for( int k = n + 2; k < n2; ++k )
    {
    output( k ) = 0;
    }

  nr::fourier::FFT< T >::four1( output.data_block(), n, -1 );

  for( int k = 0; k < n2; ++k )
    {
    output( k ) /= static_cast< T >( n );
    }
}

template< class T >
void Phase( const vnl_vector< T >& input, vnl_vector< T >& output )
{
  const int n = input.size() / 2;

  output.set_size( n );
  for( int i = 0; i < n; ++i )
    {
    output( i ) = vcl_atan2( - input[ i * 2 + 1 ], input[ i * 2 ] );
    }
}

/**
 * Instantaneous phase of a complex time varying signal.
 * @param input Vector of length n*2 (real/imaginary samples interleaved)
 * @param output Phase vector of length n
 * @param n Number of complex samples to calculate phase information of
 */
template< class T >
void Phase( const vnl_vector< T >& input, vnl_vector< T >& output, const int& n )
{
  output.set_size( n );
  for( int i = 0; i < n; ++i )
    {
    output( i ) = vcl_atan2( - input[ i * 2 + 1 ], input[ i * 2 ] );
    }
}

/**
 * Phase lag index, measures synchronization properties between two vectors
 * @param phase1 Instantaneous phases of vector 1
 * @param phase2 Instantaneous phases of vector 2
 * @return PLI
 */
template< class T >
T PhaseLagIndex( const vnl_vector< T >& phase1, const vnl_vector< T >& phase2 )
{
  const int n = phase1.size();

  int count1 = 0;
  int count2 = 0;
  for( int i = 0; i < n; ++i )
    {
    const T difference = phase1( i ) - phase2( i );

    if ( difference > 0 )
      {
      ++count1;
      ++count2;
      }
    else if ( difference < 0 )
      {
      ++count2;
      }
    }

  if ( count2 == 0 )
    {
    return 0;
    }

  return 2. * vnl_math_abs( 0.5 - static_cast< T >( count1 ) / static_cast< T >( count2 ) );
}

/**
 * Phase coherence measures phase synchronization between two vectors
 * @param phase1 Instantaneous phases of vector 1
 * @param phase2 Instantaneous phases of vector 2
 */
template< class T >
T PhaseCoherence( const vnl_vector< T >& phase1, const vnl_vector< T >& phase2 )
{
  const int n = phase1.size();

  T sumReal = 0;
  T sumImag = 0;

  for( int i = 0; i < n; ++i )
    {
    const T difference = phase1( i ) - phase2( i );
    sumReal += cos( difference );
    sumImag += sin( difference );
    }

  sumReal /= static_cast< T >( n );
  sumImag /= static_cast< T >( n );

  return sqrt( sumReal * sumReal + sumImag * sumImag );
}

template< class T >
T AbsoluteCoherency( const vnl_vector< T >& vector1, const vnl_vector< T >& vector2, const vnl_vector< T >& phase1, const vnl_vector< T >& phase2 )
{
  const int n = vector1.size();

  T sumReal = 0;
  T sumImag = 0;

  T sumAmp1 = 0;
  T sumAmp2 = 0;

  for( int i = 0; i < n; ++i )
    {
    const T difference = phase1( i ) - phase2( i );
    const T& A1 = vector1( i );
    const T& A2 = vector2( i );

    sumReal += A1 * A2 * cos( difference );
    sumImag += A1 * A2 * sin( difference );

    sumAmp1 += A1 * A1;
    sumAmp2 += A2 * A2;
    }

  sumReal /= static_cast< T >( n );
  sumImag /= static_cast< T >( n );
  sumAmp1 /= static_cast< T >( n );
  sumAmp2 /= static_cast< T >( n );

  T normalization = sqrt( sumAmp1 * sumAmp2 );

  return sqrt( ( sumReal / normalization ) * ( sumReal / normalization ) + ( sumImag / normalization ) * ( sumImag / normalization ) );
}

template< class T >
T ImaginaryCoherency( const vnl_vector< T >& vector1, const vnl_vector< T >& vector2, const vnl_vector< T >& phase1, const vnl_vector< T >& phase2 )
{
  const int n = vector1.size();

  T sumImag = 0;

  T sumAmp1 = 0;
  T sumAmp2 = 0;

  for( int i = 0; i < n; ++i )
    {
    const T difference = phase1( i ) - phase2( i );
    const T& A1 = vector1( i );
    const T& A2 = vector2( i );

    sumImag += A1 * A2 * sin( difference );

    sumAmp1 += A1 * A1;
    sumAmp2 += A2 * A2;
    }

  sumImag /= static_cast< T >( n );
  sumAmp1 /= static_cast< T >( n );
  sumAmp2 /= static_cast< T >( n );

  return sumImag / sqrt( sumAmp1 * sumAmp2 );
}


template< class T >
T PLI( const vnl_vector< T >& a, const vnl_vector< T >& b, const int& n )
{
  vnl_vector< T > vectorA, vectorB;
  vnl_vector< T > hilbertA, hilbertB;
  vnl_vector< T > phaseA, phaseB;

  ::nr::utility::RealToComplexVector< T >( a, vectorA, n );
  ::nr::utility::RealToComplexVector< T >( b, vectorB, n );

  Hilbert< T >( vectorA, n, hilbertA );
  Hilbert< T >( vectorB, n, hilbertB );

  Phase< T >( hilbertA, phaseA, a.size() );
  Phase< T >( hilbertB, phaseB, b.size() );

  return PhaseLagIndex< T >( phaseA, phaseB );
}

template< class T >
T PhaseCoherence( const vnl_vector< T >& a, const vnl_vector< T >& b, const int& n )
{
  vnl_vector< T > vectorA, vectorB;
  vnl_vector< T > hilbertA, hilbertB;
  vnl_vector< T > phaseA, phaseB;

  ::nr::utility::RealToComplexVector< T >( a, vectorA, n );
  ::nr::utility::RealToComplexVector< T >( b, vectorB, n );

  Hilbert< T >( vectorA, n, hilbertA );
  Hilbert< T >( vectorB, n, hilbertB );

  Phase< T >( hilbertA, phaseA, a.size() );
  Phase< T >( hilbertB, phaseB, b.size() );

  return PhaseCoherence< T >( phaseA, phaseB );
}

template< class T >
T AbsoluteCoherency( const vnl_vector< T >& a, const vnl_vector< T >& b, const int& n )
{
  vnl_vector< T > vectorA, vectorB;
  vnl_vector< T > hilbertA, hilbertB;
  vnl_vector< T > phaseA, phaseB;

  ::nr::utility::RealToComplexVector< T >( a, vectorA, n );
  ::nr::utility::RealToComplexVector< T >( b, vectorB, n );

  Hilbert< T >( vectorA, n, hilbertA );
  Hilbert< T >( vectorB, n, hilbertB );

  Phase< T >( hilbertA, phaseA, a.size() );
  Phase< T >( hilbertB, phaseB, b.size() );

  return AbsoluteCoherency< T >( a, b, phaseA, phaseB );
}

template< class T >
T ImaginaryCoherency( const vnl_vector< T >& a, const vnl_vector< T >& b, const int& n )
{
  vnl_vector< T > vectorA, vectorB;
  vnl_vector< T > hilbertA, hilbertB;
  vnl_vector< T > phaseA, phaseB;

  ::nr::utility::RealToComplexVector< T >( a, vectorA, n );
  ::nr::utility::RealToComplexVector< T >( b, vectorB, n );

  Hilbert< T >( vectorA, n, hilbertA );
  Hilbert< T >( vectorB, n, hilbertB );

  Phase< T >( hilbertA, phaseA, a.size() );
  Phase< T >( hilbertB, phaseB, b.size() );

  return ImaginaryCoherency< T >( a, b, phaseA, phaseB );
}

} // end namespace correlation

} // end namespace nr

#endif /*__nrCorrelation_h__*/
