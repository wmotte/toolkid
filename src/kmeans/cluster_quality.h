#ifndef CLUSTERIX_ALGORITHM_CLUSTER_QUALITY_H
#define CLUSTERIX_ALGORITHM_CLUSTER_QUALITY_H

#include <cmath>
#include <iterator>
#include <map>
#include <set>
#include <vector>

/// An algorithm to calculate various quality measures w.r.t. a given clustering
class ClusterQualityEvaluation: public MatrixAlgorithm
{
public:
	/// Constructs a clustering object for a given matrix pointer
	ClusterQualityEvaluation( Matrix* data = 0 ) : MatrixAlgorithm( data )
	{
	}
	/// Constructs a clustering object for a given matrix
	ClusterQualityEvaluation( Matrix& data ) : MatrixAlgorithm( data )
	{
	}

	/// Executes the algorithm on the given matrix
	void run( ResultMap& _result );
};

/// Calculates the mass fraction captured by a given clustering
/**
 * \tparam Iterator  the type of the iterator. It must be a random access iterator.
 * \param  begin   iterator to the beginning of a vector storing the cluster indices.
 *                 They must all be greater than zero.
 * \param  end     iterator to the end of the vector
 * \param  matrix  matrix holding the weights. It must be symmetric and it shouldn't
 *                 have any elements in the diagonal.
 */
template< typename Iterator >
double massFraction( Iterator begin, Iterator end, const Matrix& matrix )
{
	std::vector< double > cluster_indices( 1 );
	Matrix::const_iterator it = matrix.begin();
	double capturedMass = 0.0, totalMass = 0.0;

	cluster_indices[0] = -1;
	cluster_indices.insert( cluster_indices.end(), begin, end );

	while ( it != matrix.end() )
	{
		if ( cluster_indices[it->first.first] == cluster_indices[it->first.second] )
			capturedMass += it->second;
		totalMass += it->second;
		++it;
	}

	return capturedMass / totalMass;
}

/// Calculates the area fraction captured by a given clustering
/**
 * \tparam Iterator  the type of the iterator. It must be a random access iterator.
 * \param  begin   iterator to the beginning of a vector storing the cluster indices.
 *                 They must all be greater than zero.
 * \param  end     iterator to the end of the vector
 * \param  matrix  matrix holding the weights. It must be symmetric and it shouldn't
 *                 have any elements in the diagonal.
 */
template< typename Iterator >
double areaFraction( Iterator begin, Iterator end )
{
typedef	typename std::iterator_traits<Iterator>::value_type value_type;
	typedef std::multiset<value_type> set_type;

	set_type mset(begin, end);
	typename set_type::const_iterator it, it2;
	double den = 0.0;
	double result = 0.0;

	it = mset.begin();
	while (it != mset.end())
	{
		den += mset.count(*it);
		it = mset.upper_bound(*it);
	}
	den *= (den-1);

	it = mset.begin();
	while (it != mset.end())
	{
		size_t count = mset.count(*it);
		it2 = mset.upper_bound(*it);
		result += (count / den) * (count - 1);
		it = it2;
	}

	return result;
}

/// Calculates the modularity of a given clustering
/**
 * \tparam Iterator  the type of the iterator. It must be a random access iterator.
 * \param  begin   iterator to the beginning of a vector storing the cluster indices.
 *                 They must all be greater than zero.
 * \param  end     iterator to the end of the vector
 * \param  matrix  matrix holding the weights. It must be symmetric and it shouldn't
 *                 have any elements in the diagonal.
 */
template< typename Iterator >
double modularity( Iterator begin, Iterator end, const Matrix& matrix )
{
	if ( begin == end )
		return 0.0;

	int numTypes = (int) *std::max_element( begin, end ) + 1;
	std::vector< double > dIn( numTypes ), dOut( numTypes );
	std::vector< double > e( numTypes );

	double result = 0, m2 = 0;
	Matrix::const_iterator it = matrix.begin();

	while ( it != matrix.end() )
	{
		int c1 = (int) *( begin + it->first.first - 1 );
		int c2 = (int) *( begin + it->first.second - 1 );
		if ( c1 == c2 )
		{
			e[c1] += it->second;
		}
		dOut[c1] += it->second;
		dIn[c2] += it->second;
		m2 += it->second;
		++it;
	}

	result = 0.0;
	for ( int i = 0; i < numTypes; i++ )
	{
		result += e[i] / m2;
		result -= ( dIn[i] / m2 ) * ( dOut[i] / m2 );
	}

	return result;
}

/// Calculates the entropy of a clustering
/**
 * \tparam Iterator  the type of the iterator.
 * \param  begin   iterator to the beginning of a vector storing the cluster indices.
 *                 They must all be greater than zero.
 * \param  end     iterator to the end of the vector
 */
template< typename Iterator >
double entropy( Iterator begin, Iterator end )
{
	typedef std::map< typename std::iterator_traits< Iterator >::value_type, int > MapType;
	MapType counts;
	typename MapType::const_iterator it;
	double result = 0.0, n = 0;

	while ( begin != end )
	{
		counts[*begin]++;
		begin++;
		n++;
	}

	for ( it = counts.begin(); it != counts.end(); it++ )
	{
		double p = it->second / n;
		result -= p * std::log( p );
	}

	return result;
}

/// Calculates the perplexity of a clustering
/**
 * \tparam Iterator  the type of the iterator.
 * \param  begin   iterator to the beginning of a vector storing the cluster indices.
 *                 They must all be greater than zero.
 * \param  end     iterator to the end of the vector
 */
template< typename Iterator >
double perplexity( Iterator begin, Iterator end )
{
	return std::exp( entropy( begin, end ) );
}

#endif

