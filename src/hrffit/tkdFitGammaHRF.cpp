#include "tkdFitGammaHRF.h"
#include "tkdGammaHRF.h"
#include "vnl/algo/vnl_amoeba.h"

namespace tkd
{

FitGammaHRF::FitGammaHRF( const ScalarType& baseline, const int& start, const int& end )
  : m_Baseline( baseline ), m_Start( start ), m_End( end )
{
}

FitGammaHRF::ScalarType FitGammaHRF::Fit( const VectorType& data, VectorType& x, VectorType& fit, std::ostream& log )
{
  typedef GammaHRF FunctionType;
  typedef vnl_amoeba OptimizerType;

  FunctionType function( data, m_Baseline, m_Start, m_End == -1 ? data.size() : m_End );
  function.SetLogStream( log );
  OptimizerType optimizer( function );
//  optimizer.set_max_iterations( 1000000000 );
  optimizer.minimize( x );

  fit.set_size( data.size() );
  for( int i = 0; i < data.size(); ++i )
    {
    fit( i ) = function.Evaluate( i - m_Start, x );
    }

  return function.f( x );
}

} // end namespace tkd

