#include "rbmPlotData.h"

namespace rbm
{

PlotData::PlotData( VectorType data, VectorType sd ) : m_Data( data ), m_SD( sd )
{
}

QwtData* PlotData::copy() const
{
	return new PlotData( m_Data, m_SD );
}

size_t PlotData::size() const
{
	return m_Data.size();
}

double PlotData::x( size_t i ) const
{
	return i;
}

double PlotData::y( size_t i ) const
{
	return m_Data[ i ];
}

double PlotData::sd( size_t i ) const
{
	return m_SD[ i ];
}

QwtDoubleRect PlotData::boundingRect() const
{
	const size_t sz = size();

	if (sz <= 0)
		return QwtDoubleRect(1.0, 1.0, -2.0, -2.0); // invalid

	double minX, maxX, minY, maxY;
	minX = maxX = x(0);
	minY = maxY = y(0);

	for (size_t i = 1; i < sz; i++)
		{
		const double xv = x(i);
		if (xv < minX)
			minX = xv;
		if (xv > maxX)
			maxX = xv;

		const double yv = y(i) + sd(i);
		if (yv < minY)
			minY = yv;
		if (yv > maxY)
			maxY = yv;

		const double zv = y(i) - sd(i);
		if (zv < minY)
			minY = zv;
		if (zv > maxY)
			maxY = zv;
		}

	return QwtDoubleRect(minX, minY, maxX - minX, maxY - minY);

//	return QwtData::boundingRect();
}

} // end namespace rbm
