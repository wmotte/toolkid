#include "ConstrainedQuadraticProgramming.h"

/**
 *
 */
namespace cqp
{
	/**
	 * Constructor.
	 */
	ConstrainedQuadraticProgramming::ConstrainedQuadraticProgramming( const MatrixType& Q, const VectorType& q, const MatrixType& A,
			const VectorType& b, const MatrixType& C, const VectorType& d, unsigned int max_iter )
	{
		m_Q = Q;
		m_q = q;
		m_A = A;
		m_b = b;
		m_C = C;
		m_d = d;
		m_max_iter = max_iter;

		AllocateData();
	}

	/**
	 * Free mem.
	 */
	void ConstrainedQuadraticProgramming::FreeMem()
	{
		gsl_matrix_free( m_data->Q );
		gsl_vector_free( m_data->q );
		gsl_matrix_free( m_data->A );
		gsl_vector_free( m_data->b );
		gsl_matrix_free( m_data->C );
		gsl_vector_free( m_data->d );

		free( m_data );
	}

	/**
	 * Solve and return x.
	 */
	VectorType ConstrainedQuadraticProgramming::Solve()
	{
		VectorType x( m_data->Q->size1, 0 );

		unsigned int iter = 1;

		gsl_cqpminimizer *s = gsl_cqpminimizer_alloc( gsl_cqpminimizer_mg_pdip, m_data->Q->size1, m_data->A->size1, m_data->C->size1 );

		int status = gsl_cqpminimizer_set( s, m_data );

		do
		{
			status = gsl_cqpminimizer_iterate( s );
			status = gsl_cqpminimizer_test_convergence( s, 1e-10, 1e-10 );

			if ( status == GSL_SUCCESS )
			{
				for ( unsigned int j = 0; j < gsl_cqpminimizer_x( s )->size; j++ )
					x( j ) = gsl_vector_get( gsl_cqpminimizer_x( s ), j );
			} else
				iter++;
		} while ( status == GSL_CONTINUE && iter <= m_max_iter );

		// CLEAN UP ...
		gsl_cqpminimizer_free( s );

		FreeMem();

		return x;
	}

	/**
	 * Allocate data.
	 */
	void ConstrainedQuadraticProgramming::AllocateData()
	{
		// ALLOCATE ...

		m_data = static_cast< gsl_cqp_data* > ( malloc( sizeof( gsl_cqp_data ) ) );

		m_data->Q = gsl_matrix_alloc( m_Q.rows(), m_Q.cols() );
		m_data->q = gsl_vector_alloc( m_q.size() );

		m_data->A = gsl_matrix_alloc( m_A.rows(), m_A.cols() );
		m_data->b = gsl_vector_alloc( m_b.size() );

		m_data->C = gsl_matrix_alloc( m_C.rows(), m_C.cols() );
		m_data->d = gsl_vector_alloc( m_d.size() );

		// FILL

		// Q
		for( unsigned int r = 0; r < m_Q.rows(); r++ )
			for( unsigned int c = 0; c < m_Q.cols(); c++ )
				gsl_matrix_set( m_data->Q, r, c, m_Q( r, c ) );

		// q
		for( unsigned int i = 0; i < m_q.size(); i++ )
			gsl_vector_set( m_data->q, i, m_q( i ) );

		// A
		for( unsigned int r = 0; r < m_A.rows(); r++ )
			for( unsigned int c = 0; c < m_A.cols(); c++ )
				gsl_matrix_set( m_data->A, r, c, m_A( r, c ) );

		// b
		for( unsigned int i = 0; i < m_b.size(); i++ )
			gsl_vector_set( m_data->b, i, m_b( i ) );

		// C
		for( unsigned int r = 0; r < m_C.rows(); r++ )
			for( unsigned int c = 0; c < m_C.cols(); c++ )
				gsl_matrix_set( m_data->C, r, c, m_C( r, c ) );

		// d
		for( unsigned int i = 0; i < m_d.size(); i++ )
			gsl_vector_set( m_data->d, i, m_d( i ) );

		/* TEST DATA

		// objective function: 0.5*(x^t)Qx+(q^t)x
		m_data->Q = gsl_matrix_alloc( 2, 2 );
		m_data->q = gsl_vector_alloc( 2 );

		// constraints: Ax=b; Cx>=d
		m_data->A = gsl_matrix_alloc( 1, 2 );
		m_data->b = gsl_vector_alloc( 1 );
		m_data->C = gsl_matrix_calloc( 2, 2 );
		m_data->d = gsl_vector_calloc( 2 );

		// Q
		gsl_matrix_set_identity( m_data->Q );
		// q
		gsl_vector_set( m_data->q, 0, -2.0 );
		gsl_vector_set( m_data->q, 1, -1.0 );
		// A
		gsl_matrix_set( m_data->A, 0, 0, 3.0 );
		gsl_matrix_set( m_data->A, 0, 1, 1.0 );
		// b
		gsl_vector_set( m_data->b, 0, 1.8 );
		// C
		gsl_matrix_set( m_data->C, 0, 0, 1.0 );
		gsl_matrix_set( m_data->C, 1, 1, 1.0 );
		// d=0
		gsl_vector_set( m_data->d, 0, 0 );
		gsl_vector_set( m_data->d, 1, 0 );

		*/

	}
} // end namespace cqp.

