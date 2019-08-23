#ifndef __tkdGammaHRF_h__
#define __tkdGammaHRF_h__
#include "vnl/vnl_vector.h"
#include "vnl/vnl_cost_function.h"
#include <iostream>

namespace tkd
{

class GammaHRF
  : public vnl_cost_function
{
public:
  typedef double ScalarType;
  typedef vnl_vector< ScalarType > VectorType;

  GammaHRF( const VectorType& y, const ScalarType& baseline = 0, const int& start = 0, const int& end = -1 );
  ScalarType TimeToPeak( const VectorType& x ) const;
  ScalarType ChiSquared( const VectorType& x ) const;
  ScalarType ChiSquaredCDF( const ScalarType& x2 ) const;
  ScalarType Evaluate( const ScalarType& t, const VectorType& x ) const;
  ScalarType AUC( const ScalarType& t, const VectorType& x ) const;
  ScalarType InverseCDF( const ScalarType& p, const VectorType& x ) const;
  virtual ScalarType f( VectorType const& x );
  void SetLogStream( std::ostream& os );

protected:
  VectorType m_Y;
  ScalarType m_Baseline;
  int m_Start;
  int m_End;
  std::ostream* m_Log;
};

} // end namespace tkd

#endif /*__tkdGammaHRF_h__*/
