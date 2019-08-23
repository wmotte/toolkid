#ifndef CONSTRAINED_QUADRATIC_PROGRAMMING_H
#define CONSTRAINED_QUADRATIC_PROGRAMMING_H

extern "C"
{
	#include "gsl_cqp.h"
}

#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>

/**
 * Constrained Quadratic programming.
 *
 * Modified from: http://www.math.uni-trier.de/~huebner/software.html. cqp-1.2.tar.gz.
 */
namespace cqp
{
	typedef vnl_vector< double > VectorType;
	typedef vnl_matrix< double > MatrixType;

	/**
	 * Constrained Quadratic programming wrapper.
	 */
	class ConstrainedQuadraticProgramming
	{

	public:

		/**
		 * 	Constructor minimize function: 0.5x'Qx + q'x with constraints: Ax = b; Cx >= d.
		 */
		ConstrainedQuadraticProgramming( const MatrixType& Q, const VectorType& q, const MatrixType& A, const VectorType& b,
				const MatrixType& C, const VectorType& d, unsigned int max_iter );

		/**
		 * Free mem.
		 */
		void FreeMem();

		/**
		 * Solve and return x.
		 */
		VectorType Solve();

	protected:

		/**
		 * Allocate C-data.
		 */
		void AllocateData();

	private:

		gsl_cqp_data *m_data;

		MatrixType m_Q;
		VectorType m_q;
		MatrixType m_A;
		VectorType m_b;
		MatrixType m_C;
		VectorType m_d;
		unsigned int m_max_iter;

	};
} // end namespace cqp.

#endif
