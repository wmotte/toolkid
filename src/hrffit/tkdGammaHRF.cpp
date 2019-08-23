#include "tkdGammaHRF.h"
#include "nrGamma.h"
#include <iostream>

namespace tkd
{

	GammaHRF::GammaHRF( const VectorType& y, const ScalarType& baseline, const int& start, const int& end ) :
		vnl_cost_function( 6 ), m_Y( y ), m_Baseline( baseline ), m_Start( start ), m_Log( 0 )
	{
		m_End = end != -1 ? end : y.size();
	}

	void GammaHRF::SetLogStream( std::ostream& os )
	{
		m_Log = &os;
	}

	GammaHRF::ScalarType GammaHRF::Evaluate( const ScalarType& t, const VectorType& x ) const
	{
		ScalarType k = x( 0 );
		ScalarType alpha1 = x( 1 );
		ScalarType beta1 = x( 2 );

		if ( alpha1 <= 0 || beta1 <= 0 || t <= 0 || k <= 0 )
		{
			return m_Baseline;
		}

		nr::Gammadist dist( alpha1, beta1 );

		return m_Baseline + dist.p( t ) * k;//nr::Gamma::gamma( t - m_Start, alpha1, beta1 ) * k;
	}

	GammaHRF::ScalarType GammaHRF::TimeToPeak( const VectorType& x ) const
	{
		return ( x( 1 ) - 1 ) / x( 2 );
	}

	GammaHRF::ScalarType GammaHRF::ChiSquared( const VectorType& x ) const
	{
		ScalarType chiSquared = 0;

		for ( int i = m_Start + 1; i < m_End; ++i )
		{
			ScalarType t = static_cast< ScalarType > ( i - m_Start );
			ScalarType value = this->Evaluate( t, x );
			ScalarType distance = m_Y( i ) - value;
			chiSquared += ( ( distance * distance ) / value );
		}

		return chiSquared;
	}

	GammaHRF::ScalarType GammaHRF::ChiSquaredCDF( const ScalarType& x2 ) const
	{
		nr::Chisqdist dist( m_End - m_Start );
		return dist.cdf( x2 );
	}

	GammaHRF::ScalarType GammaHRF::AUC( const ScalarType& t, const VectorType& x ) const
	{
		ScalarType k = x( 0 );
		ScalarType alpha1 = x( 1 );
		ScalarType beta1 = x( 2 );

		if ( alpha1 <= 0 || beta1 <= 0 || t <= 0 || k <= 0 )
		{
			return 0;
		}

		nr::Gammadist dist( alpha1, beta1 );

		return dist.cdf( t ) * k;
	}

	GammaHRF::ScalarType GammaHRF::InverseCDF( const ScalarType& p, const VectorType& x ) const
	{
		ScalarType k = x( 0 );
		ScalarType alpha1 = x( 1 );
		ScalarType beta1 = x( 2 );

		if ( alpha1 <= 0 || beta1 <= 0 || k <= 0 )
		{
			return 0;
		}

		nr::Gammadist dist( alpha1, 1. / beta1 );

		return dist.invcdf( p );
	}

	GammaHRF::ScalarType GammaHRF::f( VectorType const& x )
	{
		ScalarType sum = 0;
		for ( int i = m_Start + 1; i < m_End; ++i )
		{
			ScalarType t = static_cast< ScalarType > ( i - m_Start );
			ScalarType value = this->Evaluate( t, x );
			ScalarType distance = value - m_Y( i );

			if ( x( 0 ) <= 0 || x( 1 ) <= 0 || x( 2 ) <= 0 )
			{
				distance = 1e2;
			}

			sum += ( distance * distance );
		}

		if ( m_Log )
		{
			( *m_Log ) << sum << "\t" << x << std::endl;
		}

		return sum;
	}

} // end namespace tkd
