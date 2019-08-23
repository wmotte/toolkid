#ifndef __dkiOptimizeFunction_h__
#define __dkiOptimizeFunction_h__

#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_cost_function.h"

#include <iostream>

namespace dki
{
	/**
	 * Diffusion kurtosis imaging functor.
	 */
	class DKIOptimizeFunction: public vnl_cost_function
	{
	public:
		typedef double ScalarType;
		typedef vnl_vector< ScalarType > VectorType;
		typedef vnl_matrix< ScalarType > MatrixType;

		/**
		 * Constructor.
		 */
		DKIOptimizeFunction( const MatrixType& A, const VectorType& b );

		/**
		 * Minimization function of x.
		 */
		virtual ScalarType f( VectorType const& x );

	protected:

		MatrixType m_A;
		VectorType m_b;
	};

} // end namespace dki

#endif /*__dkiOptimizeFunction_h__*/
