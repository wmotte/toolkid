#include "powerlawCommon.h"

/**
 * @author: W.M. Otte (wim@invivonmr.uu.nl); Image Sciences Institute, UMC Utrecht, NL.
 * @date: 19-11-2009
 *
 * Function implementations of powerlaw scaling parameter estimation.
 *
 * ***************************************************************************
 * Method: "Power-law distributions in empirical data", Clauset et al, 2009
 * http://www.santafe.edu/~aaronc/powerlaws/
 * ***************************************************************************
 */
namespace graph
{
	/**
	 * Print message of finite size bias only once when bootstrap functions are called.
	 */
	bool finiteSizeMessagePrinted = false;

	/**
	 *
	 */
	template< class ValueType >
	void Powerlaw< ValueType >::SingleFit( const VectorType& V, VectorType& results, bool nosmall, bool finite, float startValue,
			float incrementValue, float endValue )
	{
		// determine mle-type (discrete or continuous) ...
		bool discrete = IsDiscrete( V );

		if ( discrete )
		{
			std::cout << "*** INFO ***: Discrete maximum likelihood estimation." << std::endl;
			MleInt( V, nosmall, finite, startValue, incrementValue, endValue, results );

		} else
		{
			std::cout << "*** INFO ***: Continuous maximum likelihood estimation." << std::endl;
			MleReal( V, nosmall, finite, results );
		}
	}

	/**
	 *
	 */
	template< class ValueType >
	void Powerlaw< ValueType >::BootstrapFit( const VectorType& V, VectorType& results, bool nosmall, bool finiteSize,
			ValueType startValue, ValueType incrementValue, ValueType endValue, unsigned int bootstrapIterations, bool verbose )
	{
		// determine mle-type (discrete or continuous) ...
		bool discrete = IsDiscrete( V );

		Bootstrap( V, nosmall, finiteSize, startValue, incrementValue, endValue, discrete, bootstrapIterations, results, verbose );
	}

	/**
	 * Return true of all vector elements are integers.
	 */
	template< class ValueType >
	bool Powerlaw< ValueType >::IsDiscrete( const VectorType& V )
	{
		typename VectorType::const_iterator new_end = std::find_if( V.begin(), V.end(), floating_point< ValueType >() );
		return ( V.end() == new_end );
	}

	/**
	 * Remove duplicated values from vector and return unique values vector.
	 */
	template< class ValueType >
	void Powerlaw< ValueType >::Unique( const VectorType& V, VectorType& results )
	{
		VectorType tmp;

		std::copy( V.begin(), V.end(), std::back_inserter( tmp ) );

		// first sort...
		std::sort( tmp.begin(), tmp.end() );

		// determine end point of unique vector elements...
		typename VectorType::iterator new_end = std::unique( tmp.begin(), tmp.end() );

		// iterator...
		typename VectorType::iterator it = tmp.begin();

		// copy all unique elements in return vector...
		while ( it != new_end )
		{
			results.push_back( *it );
			it++;
		}
	}

	/**
	 * Remove last element of given vector.
	 */
	template< class ValueType >
	void Powerlaw< ValueType >::RemoveLastElement( VectorType& V )
	{
		V.erase( V.end() - 1, V.end() );
	}

	/**
	 * Return sorted inport vector.
	 */
	template< class ValueType >
	void Powerlaw< ValueType >::Sort( const VectorType& V, VectorType& W )
	{
		std::copy( V.begin(), V.end(), std::back_inserter( W ) );
		std::sort( W.begin(), W.end() );
	}

	/**
	 *
	 */
	template< class ValueType >
	void Powerlaw< ValueType >::KeepHigherOrEqual( VectorType& V, ValueType x )
	{
		typename VectorType::iterator new_end = std::remove_if( V.begin(), V.end(), std::bind2nd( std::less< ValueType >(), x ) );

		V.erase( new_end, V.end() );
	}

	/**
	 *
	 */
	template< class ValueType >
	void Powerlaw< ValueType >::KeepLowerOrEqual( VectorType& V, ValueType x )
	{
		typename VectorType::iterator new_end = std::remove_if( V.begin(), V.end(), std::bind2nd( std::greater_equal< ValueType >(), x ) );

		V.erase( new_end, V.end() );
	}

	/**
	 * Return incremental vector.
	 */
	template< class ValueType >
	void Powerlaw< ValueType >::GetIncrementVector( ValueType start, ValueType increment, ValueType end, VectorType& V )
	{
		int n = ( end - start ) / increment;

		if ( n <= 0 )
		{
			std::cerr << "*** WARNING ***: Increment vector is set to size: 0!" << std::endl;
		}

		for ( ; start <= end; start += increment )
		{
			V.push_back( start );
		}
	}

	/**
	 * Return standard deviation (sqrt variance).
	 */
	template< class ValueType >
	ValueType Powerlaw< ValueType >::GetSD( const std::vector< ValueType >& V )
	{
		using namespace boost;
		using namespace accumulators;

		accumulator_set< ValueType, stats< tag::variance > > acc;

		for ( unsigned int i = 0; i < V.size(); i++ )
			acc( V[i] );

		return sqrt( variance( acc ) );
	}

	/**
	 * Return cumulative sum of input vector.
	 */
	template< class ValueType >
	void Powerlaw< ValueType >::CumulativeSum( const VectorType& V, VectorType& W )
	{
		std::copy( V.begin(), V.end(), std::back_inserter( W ) );

		for ( unsigned int i = 1; i < V.size(); i++ )
		{
			W[i] += W[i - 1];
		}
	}

	/**
	 * Return uniform random values vector from inputs, with similar size.
	 */
	template< class ValueType >
	void Powerlaw< ValueType >::GetRandomValue( const VectorType& inputs, random_number_type& generator, VectorType& results )
	{
		// uniform distribution: 0 -> largest index ...
		int_distribution_type int_uni_dist( 0, inputs.size() - 1 );

		// generator ...
		int_generator_type int_distribution( generator, int_uni_dist );

		for ( unsigned int i = 0; i < inputs.size(); i++ )
		{
			results.push_back( inputs[int_distribution()] );
		}
	}

	/**
	 * Discrete Maximum likelihood estimation.
	 */
	template< class ValueType >
	void Powerlaw< ValueType >::MleInt( const VectorType& x, bool nosmall, bool finiteSize, ValueType startValue, ValueType increment,
			ValueType endValue, VectorType& results )
	{
		VectorType vec;
		GetIncrementVector( startValue, increment, endValue, vec );

		VectorType zvec( vec.size() );
		std::transform( vec.begin(), vec.end(), zvec.begin(), zeta< ValueType > () );

		VectorType xmins;
		Unique( x, xmins );

		RemoveLastElement( xmins );

		// first and second column of data matrix...
		VectorType dat1( xmins.size() );
		VectorType dat2( xmins.size() );

		VectorType sorted_x;

		Sort( x, sorted_x );

		ValueType Y = 0;
		ValueType xmax = *( std::max_element( sorted_x.begin(), sorted_x.end() ) );

		for ( unsigned int xm = 0; xm < xmins.size(); xm++ )
		{
			ValueType xmin = xmins[xm];

			VectorType z( sorted_x );

			KeepHigherOrEqual( z, xmin );

			ValueType n = z.size();

			// fill L with -Inf
			VectorType L( vec.size(), -std::numeric_limits< ValueType >::infinity() );

			// use copy of z, because z is used again later...
			VectorType tmp( z );
			std::transform( z.begin(), z.end(), tmp.begin(), log< ValueType > () );
			ValueType slogz = std::accumulate( tmp.begin(), tmp.end(), static_cast< ValueType > ( 0 ) );

			// xminvec = (1:xmin-1) ...
			VectorType xminvec_root( xmin - 1 );
			boost::iota( xminvec_root.begin(), xminvec_root.end(), 1 ); // FIXME Problem with large numbers!!

			//std::cout << "xminvec_root: ";
			//for( unsigned int i = 0; i < xminvec_root.size(); i++ )
			//	std::cout << xminvec_root.at( i ) << " ";
			//std::cout << std::endl;


			for ( unsigned int k = 0; k < vec.size(); k++ )
			{
				VectorType xminvec( xminvec_root );

				ValueType exp = vec[k];

				// xminvec.^-vec(k)...
				std::transform( xminvec.begin(), xminvec.end(), xminvec.begin(), std::bind2nd( power< ValueType > (), -exp ) );

				// sum( xminvec.^-vec(k))...
				ValueType sum = std::accumulate( xminvec.begin(), xminvec.end(), static_cast< ValueType > ( 0 ) );

				// log-likelihood
				L[k] = -exp * slogz - n * std::log( zvec[k] - sum );
			}

			typename VectorType::iterator max_it = std::max_element( L.begin(), L.end() );

			Y = *( max_it );
			unsigned int I = std::distance( L.begin(), max_it );

			//  compute KS statistic

			ValueType exp = vec[I];

			// xmin:xmax ...
			VectorType xmin_xmax( xmax - ( xmin - 1 ) );
			boost::iota( xmin_xmax.begin(), xmin_xmax.end(), xmin );

			// first_part = ( ( ( xmin:xmax ).^-exp ) ) ...
			VectorType first_part( xmin_xmax.size() );

			std::transform( xmin_xmax.begin(), xmin_xmax.end(), first_part.begin(), std::bind2nd( power< ValueType > (), -exp ) );

			// second_part = (zvec(I) - sum( ( 1:xmin-1 ).^-exp ) ) ...

			// 1:xmin - 1 ...
			VectorType pp( xmin );
			boost::iota( pp.begin(), pp.end(), 1 );
			RemoveLastElement( pp );

			// ( 1:xmin - 1 ) .^ - exp ...
			std::transform( pp.begin(), pp.end(), pp.begin(), std::bind2nd( power< ValueType > (), -exp ) );

			// sum( ( 1:xmin -1 ) .^ - exp ) ...
			ValueType sum_sp = std::accumulate( pp.begin(), pp.end(), static_cast< ValueType > ( 0 ) );

			// zvec[I] - sum( ( 1:xmin -1 ) .^ - exp ) ...
			ValueType second_part = zvec[I] - sum_sp;

			// first_part / second_part ...
			std::transform( first_part.begin(), first_part.end(), first_part.begin(), std::bind2nd( std::divides< ValueType >(),
					second_part ) );

			// fit = cumsum( first_part /. second_part ); ==
			// fit = cumsum((((xmin:xmax).^-vec(I)))./ (zvec(I) - sum((1:xmin-1).^-vec(I)))) ...
			VectorType fit;
			CumulativeSum( first_part, fit );

			// normalized histogram ...
			Histogram< ValueType > hist( z, xmin, xmax, ( xmax - xmin ), true );
			VectorType histogram = hist.getHistogram();

			//  cdi = cumsum(hist(z, xmin:xmax)./n) ...
			VectorType cdi;
			CumulativeSum( histogram, cdi );

			std::transform( fit.begin(), fit.end(), cdi.begin(), fit.begin(), std::minus< ValueType >() );
			std::transform( fit.begin(), fit.end(), fit.begin(), abs< ValueType > () );

			// dat(xm,:) = [max(abs( fit - cdi )) vec(I)] ...
			dat1[xm] = *( std::max_element( fit.begin(), fit.end() ) );
			dat2[xm] = vec[I];
		}

		std::cout << "Reached 335..." << std::endl;

		// select the index for the minimum value of D
		// [ D, I ] = min( dat( :, 1 ) ) ...
		typename VectorType::iterator min_it = std::min_element( dat1.begin(), dat1.end() );

		unsigned int I = std::distance( dat1.begin(), min_it );
		ValueType xmin = xmins[I];

		// z = x(x>=xmin);
		// n = length(z);
		VectorType z( x );
		KeepHigherOrEqual( z, xmin );
		ValueType n = z.size();

		// alpha = dat( I, 2 ) ...
		ValueType alpha = dat2[I];

		// finite-size correction
		if ( finiteSize )
		{
			alpha = alpha * ( static_cast< ValueType > ( n ) - 1 ) / static_cast< ValueType > ( n ) + 1 / static_cast< ValueType > ( n );
		}

		if ( !finiteSize && ( n < 50 ) && ( !finiteSizeMessagePrinted ) )
		{
			std::cout << "*** WARNING ***: finite-size bias may be present!" << std::endl;
			finiteSizeMessagePrinted = true;
		}

		// L = -alpha * sum( log( z ) ) - n * log( zvec( find( vec <= alpha, 1 , 'last' ) ) - sum((1:xmin-1).^-alpha));
								
		// 1:xmin - 1 ...
		VectorType pp( xmin );
		boost::iota( pp.begin(), pp.end(), 1 );
		RemoveLastElement( pp );
		// ( 1:xmin - 1 ) .^ - alpha ...
		std::transform( pp.begin(), pp.end(), pp.begin(), std::bind2nd( power< ValueType >(), -alpha ) );
		// sum( ( 1:xmin -1 ) .^ - alpha ) ...
		ValueType sum_third_part = std::accumulate( pp.begin(), pp.end(), static_cast< ValueType > ( 0 ) );
		
		std::reverse( vec.begin(), vec.end() );  
		typename VectorType::iterator it_last = std::find_if( vec.begin(), vec.end(), std::bind2nd( std::less_equal< ValueType >(), alpha ) );
		unsigned int index = std::distance( vec.begin(), it_last ) + 1;
		// n * log( zvec( find( vec <= alpha, 1 , 'last' ) ) - sum((1:xmin-1).^-alpha)) ...
		ValueType sum_second_part = n * std::log( zvec[ vec.size() - index ] - sum_third_part );
				
		// -alpha * sum( log( z ) ) ...
		std::transform( z.begin(), z.end(), z.begin(), log< ValueType >() );
		ValueType sum_first_part = -alpha * std::accumulate( z.begin(), z.end(), static_cast< ValueType > ( 0 ) );
		
		results.push_back( alpha );
		results.push_back( xmin );
		results.push_back( sum_first_part - sum_second_part );
	}

	/**
	 * Continues Maximum likelihood estimation.
	 */
	template< class ValueType >
	void Powerlaw< ValueType >::MleReal( const VectorType& x, bool nosmall, bool finiteSize, VectorType& results )
	{
		VectorType xmins;
		Unique( x, xmins );

		RemoveLastElement( xmins );

		VectorType dat( xmins.size(), 0 );

		VectorType sorted_x;

		Sort( x, sorted_x );

		for ( unsigned int xm = 0; xm < xmins.size(); xm++ )
		{
			VectorType z( sorted_x );

			ValueType xmin = xmins[xm];

			KeepHigherOrEqual( z, xmin );
			VectorType tmp( z ); // backup for later in fuction...

			ValueType n = z.size();

			// estimate alpha using direct MLE
			std::transform( z.begin(), z.end(), z.begin(), std::bind2nd( log_div< ValueType > (), xmin ) );

			ValueType a = static_cast< ValueType > ( n ) / std::accumulate( z.begin(), z.end(), static_cast< ValueType > ( 0 ) );

			if ( nosmall )
			{
				if ( ( static_cast< ValueType > ( a ) - 1 ) / sqrt( static_cast< ValueType > ( n ) ) > 0.1 )
				{
					dat.erase( dat.begin() + xm, dat.end() );
					xm = xmins.size() + 1;
					break;
				}
			}

			// compute KS statistic
			VectorType cx( n );
			boost::iota( cx.begin(), cx.end(), 0 );
			std::transform( cx.begin(), cx.end(), cx.begin(), std::bind2nd( std::divides< ValueType >(), n ) );

			// cf = xmin / z ...
			VectorType cf( tmp.size(), xmin );
			std::transform( cf.begin(), cf.end(), tmp.begin(), cf.begin(), std::divides< ValueType >() );

			// cf = 1 - ( xmin / z ) ^ a ...
			std::transform( cf.begin(), cf.end(), cf.begin(), std::bind2nd( power_minus_one< ValueType > (), a ) );

			// max( abs( cf - cx ) ) ...
			std::transform( cf.begin(), cf.end(), cx.begin(), cf.begin(), std::minus< ValueType >() );
			dat[xm] = *( std::max_element( cf.begin(), cf.end() ) );
		}

		// D = min(dat) ...
		ValueType D = *( std::min_element( dat.begin(), dat.end() ) );

		// xmin  = xmins( find( dat <= D, 1, 'first' ) ) ...
		typename VectorType::iterator new_end = std::find_if( dat.begin(), dat.end(), std::bind2nd( std::less_equal< ValueType >(), D ) );

		ValueType xmin = xmins[std::distance( dat.begin(), new_end )];

		// z = x( x >= xmin ) ...
		VectorType z( x );
		KeepHigherOrEqual( z, xmin );

		unsigned int n = z.size();

		// alpha = 1 + n ./ sum( log(z./xmin) ) ...
		VectorType tmp( z.size(), xmin );
		std::transform( z.begin(), z.end(), tmp.begin(), z.begin(), log_div< ValueType > () );
		ValueType sum = std::accumulate( z.begin(), z.end(), static_cast< ValueType > ( 0 ) );
		ValueType alpha = 1.0 + z.size() / sum;

		// finite-size correction
		if ( finiteSize )
		{
			alpha = alpha * ( static_cast< ValueType > ( n ) - 1 ) / static_cast< ValueType > ( n ) + 1 / static_cast< ValueType > ( n );
		}

		if ( !finiteSize && ( n < 50 ) && ( !finiteSizeMessagePrinted ) )
		{
			std::cout << "*** WARNING ***: finite-size bias may be present!" << std::endl;
			finiteSizeMessagePrinted = true;
		}

		// log-likelihood: L = n*log((alpha-1)/xmin) - alpha.*sum(log(z./xmin));
		ValueType L = n * std::log( ( alpha - 1 ) / xmin ) - alpha * sum;

		results.push_back( alpha );
		results.push_back( xmin );
		results.push_back( L );
	}

	/**
	 * Run bootstrapping of mle n-times and return Vector with indices:
	 *
	 * 0: average alpha
	 * 1: average xmin
	 * 2: average L
	 *
	 * 3: sd alpha
	 * 4: sd xmin
	 * 5: sd L
	 *
	 * Empty VectorType is returned when data is corrupt.
	 */
	template< class ValueType >
	void Powerlaw< ValueType >::Bootstrap( const VectorType& inputs, bool nosmall, bool finiteSize, ValueType startValue,
			ValueType increment, ValueType endValue, bool discrete, unsigned int n, VectorType& results, bool verbose )
	{
		VectorType all_alpha;
		VectorType all_xmin;
		VectorType all_L;

		// random number ...
		random_number_type generator( time( 0 ) );

		if ( verbose )
			std::cout << "*** INFO ***: boostrapping done (%): " << std::endl;

		unsigned int successfulBootstraps = 0;

		for ( unsigned int i = 0; i < n; i++ )
		{
			if ( verbose )
			{
				std::cout << ( static_cast< ValueType > ( i ) / static_cast< ValueType > ( n ) ) * 100;
				std::cout.flush();
				std::cout << '\r';
			}

			VectorType random_inputs;
			GetRandomValue( inputs, generator, random_inputs );

			VectorType run;

			Mle( random_inputs, nosmall, finiteSize, startValue, increment, endValue, discrete, run );

			if ( !run.empty() )
			{
				all_alpha.push_back( run[0] );
				all_xmin.push_back( run[1] );
				all_L.push_back( run[2] );
				successfulBootstraps++;
			}
		}

		if ( successfulBootstraps != n )
		{
			ValueType p = ( static_cast< ValueType > ( successfulBootstraps ) / static_cast< ValueType > ( n ) ) * 100;
			std::cerr << "*** WARNING ***: bootstrapping only ran partially "
				"-> (" << p << " %)." << std::endl;
		}

		if ( !all_alpha.empty() )
		{
			ValueType average_alpha = std::accumulate( all_alpha.begin(), all_alpha.end(), static_cast< ValueType > ( 0 ) )
					/ all_alpha.size();
			ValueType average_xmin = std::accumulate( all_xmin.begin(), all_xmin.end(), static_cast< ValueType > ( 0 ) ) / all_xmin.size();
			;
			ValueType average_L = std::accumulate( all_L.begin(), all_L.end(), static_cast< ValueType > ( 0 ) ) / all_L.size();
			;

			ValueType sd_alpha = GetSD( all_alpha );
			ValueType sd_xmin = GetSD( all_xmin );
			ValueType sd_L = GetSD( all_L );

			results.push_back( average_alpha );
			results.push_back( average_xmin );
			results.push_back( average_L );

			results.push_back( sd_alpha );
			results.push_back( sd_xmin );
			results.push_back( sd_L );
		}
	}

	/**
	 * Maximum likelihood estimation.
	 *
	 * If input is not correct an empty VectorType will be returned!
	 */
	template< class ValueType >
	void Powerlaw< ValueType >::Mle( const VectorType& inputs, bool nosmall, bool finiteSize, ValueType startValue, ValueType increment,
			ValueType endValue, bool discrete, VectorType& results )
	{
		// check if at all inputs are not identical ...
		if ( std::max_element( inputs.begin(), inputs.end() ) == std::min_element( inputs.begin(), inputs.end() ) )
		{
			return;
		}

		if ( discrete )
		{
			if ( startValue > 1.0 )
				MleInt( inputs, nosmall, finiteSize, startValue, increment, endValue, results );

			else
				std::cerr << "*** ERROR ***: start-value should be higher than 1.0!" << std::endl;

		} else
			MleReal( inputs, nosmall, finiteSize, results );
	}
} // end namespace graph

template class graph::Powerlaw< float >;
template class graph::Powerlaw< double >;
