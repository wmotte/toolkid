#ifndef __rbmPlotData_h__
#define __rbmPlotData_h__
#include <qwt_data.h>
#include "vnl/vnl_vector.h"
namespace rbm
{

class PlotData : public QwtData
{
public:
  typedef vnl_vector< float > VectorType;

  PlotData( VectorType data, VectorType sd );
  virtual QwtData *copy() const;
  virtual size_t size() const;
  virtual double x( size_t i ) const;
  virtual double y( size_t i ) const;
  virtual double sd( size_t i ) const;
  virtual QwtDoubleRect boundingRect() const;

protected:
  VectorType m_Data;
  VectorType m_SD;
};

} // end namespace rbm

#endif /*__rbmPlotData_h__*/
