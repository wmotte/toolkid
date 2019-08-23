#ifndef __rbmPlotter_h__
#define __rbmPlotter_h__
#include <qwt_plot.h>
#include "rbmPlotData.h"

namespace rbm
{

class PlotCurve;

class Plotter : public QwtPlot
{
public:
	typedef PlotData::VectorType VectorType;

	Plotter( VectorType mean, VectorType sd );
	void set( VectorType mean, VectorType sd );

	PlotCurve* plot;
};

} // end namespace rbm

#endif /*__rbmPlotter_h__*/
