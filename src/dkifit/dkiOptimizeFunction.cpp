#include "dkiOptimizeFunction.h"

namespace dki
{

	/**
	 * Constructor.
	 *
	 * A = DesignMatrix
	 * b = data
	 */
	DKIOptimizeFunction::DKIOptimizeFunction( const MatrixType& A, const VectorType& b )
	: vnl_cost_function( 21 )
	{
        m_A = A;
        m_b = b;
	}

	/**
	 * DKI minimization function.
	 */
	DKIOptimizeFunction::ScalarType DKIOptimizeFunction::f( VectorType const& x )
	{
		VectorType difference = m_b - m_A * x; // TODO: output idential to ULLS!


		return dot_product( difference, difference );
		//return inner_product( difference, difference );
		//return difference.sum();
	}

} // end namespace dki
