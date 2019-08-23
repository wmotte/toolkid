#ifndef __tkdFitGammaHRF_h__
#define __tkdFitGammaHRF_h__
#include "vnl/vnl_vector.h"
#include <iostream>

namespace tkd
{

class FitGammaHRF
{
public:
  typedef double ScalarType;
  typedef vnl_vector< ScalarType > VectorType;

  FitGammaHRF( const ScalarType& baseline = 0, const int& start = 0, const int& end = -1 );
  ScalarType Fit( const VectorType& data, VectorType& x, VectorType& fit, std::ostream& log = std::cout );

protected:
  ScalarType m_Baseline;
  int m_Start;
  int m_End;
};

} // end namespace tkd

#endif /*__tkdFitGammaHRF_h__*/
