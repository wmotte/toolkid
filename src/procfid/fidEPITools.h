#ifndef __fidEPITools_h__
#define __fidEPITools_h__
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_math.h"

namespace fid
{

class EPITools
{
public:
	typedef float PrecisionType;
	typedef std::complex< PrecisionType > ComplexType;
	typedef vnl_vector< PrecisionType > VectorType;
	typedef vnl_matrix< PrecisionType > MatrixType;
	typedef vnl_vector< ComplexType > ComplexVectorType;
	typedef vnl_matrix< ComplexType > ComplexMatrixType;
	typedef vnl_vector< int > MaskVectorType;

	static VectorType CalculateWeights( const ComplexMatrixType& matrix );
	static void CalculatePhase( const ComplexVectorType& trace, VectorType& output );
	static void PhaseWithMap( ComplexVectorType& input, const VectorType& phasemap );
	static MaskVectorType BuildMask( const ComplexVectorType& trace, PrecisionType range );
	static void FitPhase( VectorType& phase, const MaskVectorType& mask, const VectorType& weights );
	static VectorType BuildPhaseMap( ComplexVectorType trace, const VectorType& weights );
};

} // end namespace fid

#endif /*__fidEPITools_h__*/
