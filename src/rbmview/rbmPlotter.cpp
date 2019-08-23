#include "rbmPlotter.h"
#include "rbmPlotData.h"
#include "rbmPlotCurve.h"
#include <qwt_legend.h>
#include <qwt_text.h>

namespace rbm
{

Plotter::Plotter( VectorType mean, VectorType sd )
{
	setTitle( "Mean time course" );

	setAxisTitle( xBottom, "volume" );
	setAxisTitle( yLeft, "intensity" );

	plot = new PlotCurve( "mean" );
	plot->attach( this );
	plot->setData( PlotData( mean, sd ) );
}

void Plotter::set( VectorType mean, VectorType sd )
{
	plot->setData( PlotData( mean, sd ) );
	replot();
}

} // end namespace rbm
