#include "fidEPITools.h"
#include "vnl/algo/vnl_svd.h"
#include <vector>
#include <algorithm>

namespace fid
{

EPITools::VectorType EPITools::CalculateWeights( const ComplexMatrixType& matrix )
{
	// calculate squared absolute values

	MatrixType A( matrix.rows(), matrix.cols() );
	for( unsigned int i = 0; i < matrix.rows(); ++i )
		{
		for( unsigned int j = 0; j < matrix.cols(); ++j )
			{
			const std::complex< float >& p = matrix( i, j );

			A( i, j ) = p.real() * p.real() + p.imag() * p.imag();
			}
		}

	// find maximum values for each column

	vnl_matrix< float > B( 1, matrix.cols() );
	B.fill( 0 );
	for( unsigned int j = 0; j < B.size(); ++j )
		{
		float max = 0;
		for( unsigned int i = 0; i < A.rows(); ++i )
			{
			max = max < A( i, j ) ? A( i, j ) : max;
			}
		B( 0, j ) = max;
		}

	// fit B' * x = A' for x

	vnl_svd< PrecisionType > svd( B.transpose() );

	return svd.solve( A.transpose() ).get_row( 0 );
}

void EPITools::CalculatePhase( const ComplexVectorType& trace, VectorType& output )
{
	VectorType phase1 = VectorType( trace.size() );
	VectorType phase2 = VectorType( trace.size() );

	for( unsigned int i = 0; i < phase1.size(); ++i )
		{
		phase1( i ) = atan2( trace( i ).real(), trace( i ).imag() );

		// phase2 = mod( phase1 + 2PI, 2PI )
		static const PrecisionType pi2 = 2. * vnl_math::pi;
		const PrecisionType a = ( phase1( i ) + pi2 );
		const PrecisionType b = a / pi2;
		phase2( i ) = a - ( vcl_floor( b ) * pi2 );

		output( i ) = phase1( i );
		}

	for( unsigned int i = 1; i < output.size(); ++i )
		{
		PrecisionType d1 = phase1( i - 1 ) - phase1( i );
		PrecisionType d2 = phase2( i - 1 ) - phase2( i );

		output( i ) = output( i - 1 ) + ( vnl_math_abs( d1 ) < vnl_math_abs( d2 ) ? d1 : d2 );
		}

	int index = round( output.size() / 2. ) - 1;
	float mean = output( index ) - phase1( index );
	for( unsigned int i = 0; i < output.size(); ++i )
		{
		output( i ) -= mean;
		}
}

void EPITools::PhaseWithMap( ComplexVectorType& input, const VectorType& phasemap )
{
	for( unsigned int i = 0; i < input.size(); ++i )
		{
		input( i ) = ComplexType(
										cos( phasemap( i ) ) * input( i ).real() + sin( phasemap( i ) ) * input( i ).imag(),
										cos( phasemap( i ) ) * input( i ).imag() - sin( phasemap( i ) ) * input( i ).real() );
		}
}

EPITools::MaskVectorType EPITools::BuildMask( const ComplexVectorType& trace, PrecisionType range )
{
	std::vector< PrecisionType > list;
	for( unsigned int i = 0; i < trace.size(); ++i )
		{
		list.push_back( vcl_sqrt( trace[ i ].real() * trace[ i ].real() + trace[ i ].imag() * trace[ i ].imag() ) );
		}

	std::vector< PrecisionType > copy = list;
	std::sort( list.begin(), list.end() );
	int indexLow = vcl_ceil( ( 1.0 - range ) * static_cast< PrecisionType >( trace.size() ) );
	int indexHigh = vcl_floor( range * static_cast< PrecisionType >( trace.size() ) );

	PrecisionType low = list[ indexLow ];
	PrecisionType high = list[ indexHigh ];

	MaskVectorType mask = MaskVectorType( trace.size() );
	for( unsigned int i = 0; i < trace.size(); ++i )
		{
		mask[ i ] = ( copy[ i ] >= low && copy[ i ] <= high ) ? 1 : 0;
		}

	return mask;
}

void EPITools::FitPhase( VectorType& phase, const MaskVectorType& mask, const VectorType& weights )
{
	int n = 0;

	for( unsigned int i = 0; i < mask.size(); ++i )
		{
		n += mask( i );
		}

	vnl_vector< double > xs = vnl_vector< double >( n );
	vnl_vector< double > ys = vnl_vector< double >( n );
	vnl_vector< double > ws = vnl_vector< double >( n );

	n = 0;
	for( unsigned int i = 0; i < mask.size(); ++i )
		{
		if ( mask( i ) == 1 )
			{
			xs( n ) = static_cast< double >( i );
			ys( n ) = phase( i );
			ws( n ) = weights( i );
			++n;
			}
		}

	int m = 3;

	// y = a * x^2 + b * x + c;
	vnl_matrix< double > A( n, m );
	vnl_vector< double > b( n );

	for( int i = 0; i < n; ++i )
		{
		A( i, 0 ) = xs( i ) * xs( i ) * ws( i );
		A( i, 1 ) = xs( i ) * ws( i );
		A( i, 2 ) = ws( i );
		ys( i ) *= ws( i );
		}

	vnl_svd< double > svd( A );
	vnl_matrix< double > U = svd.U();
	vnl_matrix< double > W = svd.W();
	vnl_matrix< double > V = svd.V();

	vnl_vector< double > a( m );
	a.fill( 0 );

	for( int j = 0; j < m; ++j )
		{
		for( int i = 0; i < m; ++i )
			{
			double dot = 0;
			for( int k = 0; k < n; ++k )
				{
				dot += U( k, i ) * ys( k );
				}

			double v = ( dot / W( i, i ) ) * V( j, i );

			a( j ) += v;
			}
		}

	for( unsigned int i = 0; i < phase.size(); ++i )
		{
		double x = i;

		phase( i ) = static_cast< PrecisionType >( a( 0 ) * x * x + a( 1 ) * x + a( 2 ) );
		}
}

EPITools::VectorType EPITools::BuildPhaseMap( ComplexVectorType trace, const VectorType& weights )
{
	unsigned int length = trace.size();

	VectorType phase( length );
	VectorType phase2( length );

	// mask
	MaskVectorType mask = BuildMask( trace, 0.95 );

	// calculate phase
	CalculatePhase( trace, phase );

	// fit polynomial
	FitPhase( phase, mask, weights );

	// map
	PhaseWithMap( trace, phase );

	// calculate phase
	CalculatePhase( trace, phase2 );

	// store
	for( unsigned int x = 0; x < length; ++x )
		{
		phase( x ) -= phase2( round( length / 2. ) - 1 );
		}

	return phase;
}

} // end namespace fid
