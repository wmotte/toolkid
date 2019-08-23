#include "fidFourier.h"
#include "itkNumericTraits.h"
#include "vnl/algo/vnl_fft_1d.h"

namespace fid
{

template< class TPrecision >
struct Fourier< TPrecision >::Impl
{
	int N;
	int M;
	bool OutputKSpace;
};

template< class TPrecision >
Fourier< TPrecision >::Fourier( int N, int M )
: m_Impl( new Impl, true )
{
	m_Impl->N = N;
	m_Impl->M = M;
	m_Impl->OutputKSpace = false;
}

template< class TPrecision >
Fourier< TPrecision >::~Fourier()
{
}

template< class TPrecision >
void Fourier< TPrecision >::SetOutputKSpace( bool k )
{
  m_Impl->OutputKSpace = k;
}

template< class TPrecision >
bool Fourier< TPrecision >::GetOutputKSpace() const
{
  return m_Impl->OutputKSpace;
}

template< class TPrecision >
void Fourier< TPrecision >::ifft( VectorType& v, int N )
{
  if ( static_cast< int >( v.size() ) > N )
    {
    double length = v.size() - N;
    double half = round( length / 2. );
    int start = half;
    int end = half + N;

    VectorType w( N );
    for( int i = start; i < end; ++i )
      {
      w( i - start ) = v( i );
      }

    v = w;
    }

	const int n = static_cast< int >( v.size() );
	const PrecisionType factor = 1. / static_cast< PrecisionType >( N );

  VectorType w( N );
  w.fill( itk::NumericTraits< PixelType >::Zero );

  const int shift1 = N - static_cast< int >( round( static_cast< float >( n ) / 2. ) );
	for( int i = 0; i < n; ++i )
		{
		const PixelType& p = v( i );
		int index = ( i + shift1 ) % N;
		w( index ) = PixelType( p.imag(), p.real() ) * factor;
		}

	if ( !m_Impl->OutputKSpace )
	  {
	  vnl_fft_1d< PrecisionType > fft( N );
	  fft.bwd_transform( w );
	  }

	const int shift2 = N / 2;
	if ( n < N )
		{
		v.set_size( N );
		}

	for( int i = 0; i < N; ++i )
		{
		v( ( i + shift2 ) % N ) = PixelType( w( i ).imag(), w( i ).real() );
		}
}

template< class TPrecision >
void Fourier< TPrecision >::ifft( VectorType& vector )
{
	return ifft( vector, m_Impl->N );
}

template< class TPrecision >
void Fourier< TPrecision >::ifft( MatrixType& matrix )
{
	const int n = static_cast< int >( matrix.rows() );
	const int m = static_cast< int >( matrix.cols() );

	const int& N = m_Impl->N;
	const int& M = m_Impl->M;

	if ( n == N && m == M )
		{
		for( int i = 0; i < N; ++i )
			{
			VectorType v = matrix.get_row( i );
			ifft( v, M );
			matrix.set_row( i, v );
			}

		for( int i = 0; i < M; ++i )
			{
			VectorType v = matrix.get_column( i );
			ifft( v, N );
			matrix.set_column( i, v );
			}
		}
	else
		{
		MatrixType tmp( n, M );
		tmp.fill( itk::NumericTraits< PixelType >::Zero );

		for( int i = 0; i < n; ++i )
			{
			VectorType v = matrix.get_row( i );
			ifft( v, M );
			tmp.set_row( i, v );
			}

		matrix.set_size( N, M );

		for( int i = 0; i < M; ++i )
			{
			VectorType v = tmp.get_column( i );
			ifft( v, N );
			matrix.set_column( i, v );
			}
		}

}

} // end namespace fid

template class fid::Fourier< float >;
