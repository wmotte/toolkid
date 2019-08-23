#ifndef __rsPCA_h__
#define __rsPCA_h__
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "vnl/algo/vnl_svd.h"
#include "vcl_cmath.h"
#include <vector>

namespace rs
{

class PCA
{
public:
  typedef float PixelType;
  typedef vnl_matrix< float > MatrixType;
  typedef vnl_vector< float > VectorType;

  PCA( const MatrixType& x );
  void TransformPCA();
  unsigned int GetNumberOfComponents() const;
  VectorType GetComponent( unsigned int i ) const;
  PixelType GetVariance( unsigned int i ) const;
  PixelType GetCumulativeVariance( unsigned int i ) const;

protected:
  MatrixType m_X;
  MatrixType m_D;
  VectorType m_Mean;
  MatrixType m_b;
  VectorType m_L;
  unsigned int m_Components;
  std::vector< PixelType > m_Variance;
  std::vector< PixelType > m_CumulativeVariance;
};

} // end namespace rs

#endif /*__rsPCA_h__*/
