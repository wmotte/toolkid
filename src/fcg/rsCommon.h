#ifndef __rsCommon_h__
#define __rsCommon_h__
#include "vnl/vnl_vector.h"

namespace rs
{

template< class PixelType >
class Common
{
public:
  typedef vnl_vector< PixelType > VectorType;

  static PixelType CorrelationCoefficient( const VectorType& a, const VectorType& b );
  static PixelType RtoZ( const PixelType& r );
  static PixelType ZtoR( const PixelType& z );
};


} // end namespace rs

#endif /*__rsCommon_h__*/
