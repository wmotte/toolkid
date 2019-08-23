/* -*- mode: C++ ; c-file-style: "stroustrup" -*- *****************************
 * Qwt Widget Library
 * Copyright (C) 1997   Josef Wilgen
 * Copyright (C) 2002   Uwe Rathmann
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the Qwt License, Version 1.0
 *****************************************************************************/
#include <iostream>
#include <qpainter.h>
#include <qpixmap.h>
#include <qbitarray.h>
#include "qwt_global.h"
#include "qwt_legend.h"
#include "qwt_legend_item.h"
#include "rbmPlotData.h"
#include "qwt_rect.h"
#include "qwt_scale_map.h"
#include "qwt_double_rect.h"
#include "qwt_math.h"
#include "qwt_painter.h"
#include "qwt_plot.h"
#include "qwt_plot_canvas.h"
#include "qwt_curve_fitter.h"
#include "qwt_symbol.h"
#include "rbmPlotCurve.h"

#if QT_VERSION >= 0x040000

#include <qevent.h>

namespace rbm
{

class PlotCurvePaintHelper: public QObject
{
public:
    PlotCurvePaintHelper(const PlotCurve *curve, int from, int to):
        _curve(curve),
        _from(from),
        _to(to)
    {
    }

    virtual bool eventFilter(QObject *, QEvent *event)
    {
        if ( event->type() == QEvent::Paint )
        {
            _curve->draw(_from, _to);
            return true;
        }
        return false;
    }
private:
    const PlotCurve *_curve;
    int _from;
    int _to;
};

#endif // QT_VERSION >= 0x040000

static int verifyRange(int size, int &i1, int &i2)
{
    if (size < 1)
        return 0;

    i1 = qwtLim(i1, 0, size-1);
    i2 = qwtLim(i2, 0, size-1);

    if ( i1 > i2 )
        qSwap(i1, i2);

    return (i2 - i1 + 1);
}

class PlotCurve::PrivateData
{
public:
    class PixelMatrix: private QBitArray
    {
    public:
        PixelMatrix(const QRect& rect):
            QBitArray(rect.width() * rect.height()),
            _rect(rect)
        {
            fill(false);
        }

        inline bool testPixel(const QPoint& pos)
        {
            if ( !_rect.contains(pos) )
                return false;

            const int idx = _rect.width() * (pos.y() - _rect.y()) +
                (pos.x() - _rect.x());

            const bool marked = testBit(idx);
            if ( !marked )
                setBit(idx, true);

            return !marked;
        }

    private:
        QRect _rect;
    };

    PrivateData():
        curveType(Yfx),
        style(PlotCurve::Lines),
        reference(0.0),
        attributes(0),
        paintAttributes(0)
    {
        pen = QPen(Qt::black, 0);
        curveFitter = new QwtSplineCurveFitter;
    }

    ~PrivateData()
    {
        delete curveFitter;
    }

    PlotCurve::CurveType curveType;
    PlotCurve::CurveStyle style;
    double reference;

    QwtSymbol sym;
    QwtCurveFitter *curveFitter;

    QPen pen;
    QBrush brush;

    int attributes;
    int paintAttributes;
};

/*!
  \brief Ctor
*/
PlotCurve::PlotCurve():
    QwtPlotItem(QwtText())
{
    init();
}

/*!
  \brief Ctor
  \param title title of the curve
*/
PlotCurve::PlotCurve(const QwtText &title):
    QwtPlotItem(title)
{
    init();
}

/*!
  \brief Ctor
  \param title title of the curve
*/
PlotCurve::PlotCurve(const QString &title):
    QwtPlotItem(QwtText(title))
{
    init();
}

//! Dtor
PlotCurve::~PlotCurve()
{
    delete d_xy;
    delete d_data;
}

/*!
  \brief Initialize data members
*/
void PlotCurve::init()
{
    setItemAttribute(QwtPlotItem::Legend);
    setItemAttribute(QwtPlotItem::AutoScale);

    d_data = new PrivateData;
    d_xy = new PlotData( PlotData::VectorType(), PlotData::VectorType() );

    setZ(20.0);
}

int PlotCurve::rtti() const
{
    return QwtPlotItem::Rtti_PlotCurve;
}

/*!
  \brief Specify an attribute how to draw the curve

  The attributes can be used to modify the drawing algorithm.

  The following attributes are defined:<dl>
  <dt>PaintFiltered</dt>
  <dd>Tries to reduce the data that has to be painted, by sorting out
      duplicates, or paintings outside the visible area. Might have a
      notable impact on curves with many close points.
      Only a couple of very basic filtering algos are implemented.</dd>
  <dt>ClipPolygons</dt>
  <dd>Clip polygons before painting them.
  </dl>

  The default is, that no paint attributes are enabled.

  \param attribute Paint attribute
  \param on On/Off
  /sa testPaintAttribute()
*/
void PlotCurve::setPaintAttribute(PaintAttribute attribute, bool on)
{
    if ( on )
        d_data->paintAttributes |= attribute;
    else
        d_data->paintAttributes &= ~attribute;
}

/*!
    \brief Return the current paint attributes
    \sa setPaintAttribute
*/
bool PlotCurve::testPaintAttribute(PaintAttribute attribute) const
{
    return (d_data->paintAttributes & attribute);
}

/*!
  \brief Set the curve's drawing style

  Valid styles are:
  <dl>
  <dt>NoCurve</dt>
  <dd>Don't draw a curve. Note: This doesn't affect the symbol. </dd>
  <dt>Lines</dt>
  <dd>Connect the points with straight lines. The lines might
      be interpolated depending on the 'Fitted' option. Curve
      fitting can be configured using setCurveFitter.</dd>
  <dt>Sticks</dt>
  <dd>Draw vertical sticks from a baseline which is defined by setBaseline().</dd>
  <dt>Steps</dt>
  <dd>Connect the points with a step function. The step function
      is drawn from the left to the right or vice versa,
      depending on the 'Inverted' option.</dd>
  <dt>Dots</dt>
  <dd>Draw dots at the locations of the data points. Note:
      This is different from a dotted line (see setPen()).</dd>
  <dt>UserCurve ...</dt>
  <dd>Styles >= UserCurve are reserved for derived
      classes of PlotCurve that overload drawCurve() with
      additional application specific curve types.</dd>
  </dl>
  \sa style()
*/
void PlotCurve::setStyle(CurveStyle style)
{
    if ( style != d_data->style )
    {
        d_data->style = style;
        itemChanged();
    }
}

/*!
    \brief Return the current style
    \sa setStyle
*/
PlotCurve::CurveStyle PlotCurve::style() const
{
    return d_data->style;
}

/*!
  \brief Assign a symbol
  \param s symbol
  \sa symbol()
*/
void PlotCurve::setSymbol(const QwtSymbol &s )
{
    d_data->sym = s;
    itemChanged();
}

/*!
    \brief Return the current symbol
    \sa setSymbol
*/
const QwtSymbol &PlotCurve::symbol() const
{
    return d_data->sym;
}

/*!
  \brief Assign a pen
  \param p New pen
  \sa pen(), brush()
*/
void PlotCurve::setPen(const QPen &p)
{
    if ( p != d_data->pen )
    {
        d_data->pen = p;
        itemChanged();
    }
}

/*!
    \brief Return the pen used to draw the lines
    \sa setPen(), brush()
*/
const QPen& PlotCurve::pen() const
{
    return d_data->pen;
}

/*!
  \brief Assign a brush.
         In case of brush.style() != QBrush::NoBrush
         and style() != PlotCurve::Sticks
         the area between the curve and the baseline will be filled.
         In case !brush.color().isValid() the area will be filled by
         pen.color(). The fill algorithm simply connects the first and the
         last curve point to the baseline. So the curve data has to be sorted
         (ascending or descending).
  \param brush New brush
  \sa brush(), setBaseline(), baseline()
*/
void PlotCurve::setBrush(const QBrush &brush)
{
    if ( brush != d_data->brush )
    {
        d_data->brush = brush;
        itemChanged();
    }
}

/*!
  \brief Return the brush used to fill the area between lines and the baseline
  \sa setBrush(), setBaseline(), baseline()
*/
const QBrush& PlotCurve::brush() const
{
    return d_data->brush;
}

/*!
  Initialize data with a pointer to PlotData.

  \param data Data
  \sa PlotData::copy()
*/
void PlotCurve::setData(const PlotData &data)
{
    delete d_xy;
    d_xy = static_cast< PlotData* >( data.copy() );
    itemChanged();
}


QwtDoubleRect PlotCurve::boundingRect() const
{
    if ( d_xy == NULL )
        return QwtDoubleRect(1.0, 1.0, -2.0, -2.0); // invalid

    return d_xy->boundingRect();
}

/*!
  \brief Draw the complete curve

  \param painter Painter
  \param xMap Maps x-values into pixel coordinates.
  \param yMap Maps y-values into pixel coordinates.

  \sa drawCurve(), drawSymbols()
*/
void PlotCurve::draw(QPainter *painter,
    const QwtScaleMap &xMap, const QwtScaleMap &yMap,
    const QRect &) const
{
    draw(painter, xMap, yMap, 0, -1);
}

/*!
  \brief Draw a set of points of a curve.

  When observing an measurement while it is running, new points have to be
  added to an existing curve. drawCurve can be used to display them avoiding
  a complete redraw of the canvas.

  Setting plot()->canvas()->setAttribute(Qt::WA_PaintOutsidePaintEvent, true);
  will result in faster painting, if the paint engine of the canvas widget
  supports this feature.

  \param from Index of the first point to be painted
  \param to Index of the last point to be painted. If to < 0 the
         curve will be painted to its last point.

  \sa drawCurve(), drawSymbols()
*/
void PlotCurve::draw(int from, int to) const
{
    if ( !plot() )
        return;

    QwtPlotCanvas *canvas = plot()->canvas();

    bool directPaint = true;

#if QT_VERSION >= 0x040000
    if ( !canvas->testAttribute(Qt::WA_WState_InPaintEvent) &&
        !canvas->testAttribute(Qt::WA_PaintOutsidePaintEvent) )
    {
        /*
          We save curve and range in helper and call repaint.
          The helper filters the Paint event, to repeat
          the PlotCurve::draw, but now from inside the paint
          event.
         */

        PlotCurvePaintHelper helper(this, from, to);
        canvas->installEventFilter(&helper);
        canvas->repaint();

        return;
    }
#endif

    const QwtScaleMap xMap = plot()->canvasMap(xAxis());
    const QwtScaleMap yMap = plot()->canvasMap(yAxis());

    if ( canvas->testPaintAttribute(QwtPlotCanvas::PaintCached) &&
        canvas->paintCache() && !canvas->paintCache()->isNull() )
    {
        QPainter cachePainter((QPixmap *)canvas->paintCache());
        cachePainter.translate(-canvas->contentsRect().x(),
            -canvas->contentsRect().y());

        draw(&cachePainter, xMap, yMap, from, to);
    }

    if ( directPaint )
    {
        QPainter painter(canvas);

        painter.setClipping(true);
        painter.setClipRect(canvas->contentsRect());

        draw(&painter, xMap, yMap, from, to);

        return;
    }

#if QT_VERSION >= 0x040000
    if ( canvas->testPaintAttribute(QwtPlotCanvas::PaintCached) &&
        canvas->paintCache() )
    {
        /*
          The cache is up to date. We flush it via repaint to the
          canvas. This works flicker free but is much ( > 10x )
          slower than direct painting.
         */

        const bool noBG = canvas->testAttribute(Qt::WA_NoBackground);
        if ( !noBG )
            canvas->setAttribute(Qt::WA_NoBackground, true);

        canvas->repaint(canvas->contentsRect());

        if ( !noBG )
            canvas->setAttribute(Qt::WA_NoBackground, false);

        return;
    }
#endif

    // Ok, we give up
    canvas->repaint(canvas->contentsRect());
}

/*!
  \brief Draw an interval of the curve
  \param painter Painter
  \param xMap maps x-values into pixel coordinates.
  \param yMap maps y-values into pixel coordinates.
  \param from index of the first point to be painted
  \param to index of the last point to be painted. If to < 0 the
         curve will be painted to its last point.

  \sa drawCurve(), draSymbols(),
*/
void PlotCurve::draw(QPainter *painter,
    const QwtScaleMap &xMap, const QwtScaleMap &yMap,
    int from, int to) const
{
    if ( !painter || dataSize() <= 0 )
        return;

    if (to < 0)
        to = dataSize() - 1;

    if ( verifyRange(dataSize(), from, to) > 0 )
    {
        painter->save();
        painter->setPen(d_data->pen);

        /*
          Qt 4.0.0 is slow when drawing lines, but it´s even
          slower when the painter has a brush. So we don't
          set the brush before we really need it.
         */

        drawCurve(painter, d_data->style, xMap, yMap, from, to);
        painter->restore();

        if (d_data->sym.style() != QwtSymbol::NoSymbol)
        {
            painter->save();
            drawSymbols(painter, d_data->sym, xMap, yMap, from, to);
            painter->restore();
        }
    }
}

/*!
  \brief Draw the line part (without symbols) of a curve interval.
  \param painter Painter
  \param style curve style, see PlotCurve::CurveStyle
  \param xMap x map
  \param yMap y map
  \param from index of the first point to be painted
  \param to index of the last point to be painted
  \sa draw(), drawDots(), drawLines(), drawSteps(), drawSticks()
*/

void PlotCurve::drawCurve(QPainter *painter, int style,
    const QwtScaleMap &xMap, const QwtScaleMap &yMap,
    int from, int to) const
{
    switch (style)
    {
        case Lines:
            if ( testCurveAttribute(Fitted) )
            {
                // we always need the complete
                // curve for fitting
                from = 0;
                to = dataSize() - 1;
            }
            drawLines(painter, xMap, yMap, from, to);
            break;
        case Sticks:
            drawSticks(painter, xMap, yMap, from, to);
            break;
        case Steps:
            drawSteps(painter, xMap, yMap, from, to);
            break;
        case Dots:
            drawDots(painter, xMap, yMap, from, to);
            break;
        case NoCurve:
        default:
            break;
    }
}

/*!
  \brief Draw lines

  If the CurveAttribute Fitted is enabled a QwtCurveFitter tries
  to interpolate/smooth the curve, before it is painted.

  \param painter Painter
  \param xMap x map
  \param yMap y map
  \param from index of the first point to be painted
  \param to index of the last point to be painted

  \sa setCurveAttribute(), setCurveFitter(), draw(),
      drawLines(), drawDots(), drawSteps(), drawSticks()
*/
void PlotCurve::drawLines(QPainter *painter,
    const QwtScaleMap &xMap, const QwtScaleMap &yMap,
    int from, int to) const
{
    int size = to - from + 1;
    if ( size <= 0 )
        return;

    QwtPolygon polyline;
		polyline.resize(size);

		painter->setPen( QPen( Qt::red ) );

		for (int i = from; i <= to; i++)
		{
			int xi = xMap.transform(x(i));
			int yi = yMap.transform(y(i));
			int sd1 = yMap.transform(y(i)-sd(i));
			int sd2 = yMap.transform(y(i)+sd(i));
			polyline.setPoint(i - from, xi, yi);

			QwtPainter::drawLine(painter,xi,sd1,xi,sd2);
		}

		painter->setPen( QPen( Qt::black ) );
		painter->setBackground( QBrush( Qt::white ) );

		QwtPainter::drawPolyline(painter, polyline);
}

/*!
  Draw sticks

  \param painter Painter
  \param xMap x map
  \param yMap y map
  \param from index of the first point to be painted
  \param to index of the last point to be painted

  \sa draw(), drawCurve(), drawDots(), drawLines(), drawSteps()
*/
void PlotCurve::drawSticks(QPainter *painter,
    const QwtScaleMap &xMap, const QwtScaleMap &yMap,
    int from, int to) const
{
    int x0 = xMap.transform(d_data->reference);
    int y0 = yMap.transform(d_data->reference);

    for (int i = from; i <= to; i++)
    {
        const int xi = xMap.transform(x(i));
        const int yi = yMap.transform(y(i));

        if (d_data->curveType == Xfy)
            QwtPainter::drawLine(painter, x0, yi, xi, yi);
        else
            QwtPainter::drawLine(painter, xi, y0, xi, yi);
    }
}

/*!
  Draw dots

  \param painter Painter
  \param xMap x map
  \param yMap y map
  \param from index of the first point to be painted
  \param to index of the last point to be painted

  \sa draw(), drawCurve(), drawSticks(), drawLines(), drawSteps()
*/
void PlotCurve::drawDots(QPainter *painter,
    const QwtScaleMap &xMap, const QwtScaleMap &yMap,
    int from, int to) const
{
    const QRect window = painter->window();
    if ( window.isEmpty() )
        return;

    const bool doFill = d_data->brush.style() != Qt::NoBrush;

    QwtPolygon polyline;
    if ( doFill )
        polyline.resize(to - from + 1);

    if ( to > from && d_data->paintAttributes & PaintFiltered )
    {
        if ( doFill )
        {
            QPoint pp( xMap.transform(x(from)), yMap.transform(y(from)) );

            QwtPainter::drawPoint(painter, pp.x(), pp.y());
            polyline.setPoint(0, pp);

            int count = 1;
            for (int i = from + 1; i <= to; i++)
            {
                const QPoint pi(xMap.transform(x(i)), yMap.transform(y(i)));
                if ( pi != pp )
                {
                    QwtPainter::drawPoint(painter, pi.x(), pi.y());

                    polyline.setPoint(count, pi);
                    count++;

                    pp = pi;
                }
            }
            if ( int(polyline.size()) != count )
                polyline.resize(count);
        }
        else
        {
            // if we don't need to fill, we can sort out
            // duplicates independent from the order

            PrivateData::PixelMatrix pixelMatrix(window);

            for (int i = from; i <= to; i++)
            {
                const QPoint p( xMap.transform(x(i)),
                    yMap.transform(y(i)) );

                if ( pixelMatrix.testPixel(p) )
                    QwtPainter::drawPoint(painter, p.x(), p.y());
            }
        }
    }
    else
    {
        for (int i = from; i <= to; i++)
        {
            const int xi = xMap.transform(x(i));
            const int yi = yMap.transform(y(i));
            QwtPainter::drawPoint(painter, xi, yi);

            if ( doFill )
                polyline.setPoint(i - from, xi, yi);
        }
    }

    if ( doFill )
    {
        if ( d_data->paintAttributes & ClipPolygons )
        {
            const QwtRect r = painter->window();
            polyline = r.clip(polyline);
        }

        fillCurve(painter, xMap, yMap, polyline);
    }
}

/*!
  Draw step function

  The direction of the steps depends on Inverted attribute.

  \param painter Painter
  \param xMap x map
  \param yMap y map
  \param from index of the first point to be painted
  \param to index of the last point to be painted

  \sa CurveAttribute, setCurveAttribute(),
      draw(), drawCurve(), drawDots(), drawLines(), drawSticks()
*/
void PlotCurve::drawSteps(QPainter *painter,
    const QwtScaleMap &xMap, const QwtScaleMap &yMap,
    int from, int to) const
{
    QwtPolygon polyline(2 * (to - from) + 1);

    bool inverted = d_data->curveType == Yfx;
    if ( d_data->attributes & Inverted )
        inverted = !inverted;

    int i,ip;
    for (i = from, ip = 0; i <= to; i++, ip += 2)
    {
        const int xi = xMap.transform(x(i));
        const int yi = yMap.transform(y(i));

        if ( ip > 0 )
        {
            if (inverted)
                polyline.setPoint(ip - 1, polyline[ip-2].x(), yi);
            else
                polyline.setPoint(ip - 1, xi, polyline[ip-2].y());
        }

        polyline.setPoint(ip, xi, yi);
    }

    if ( d_data->paintAttributes & ClipPolygons )
    {
        const QwtRect r = painter->window();
        polyline = r.clip(polyline);
    }

    QwtPainter::drawPolyline(painter, polyline);

    if ( d_data->brush.style() != Qt::NoBrush )
        fillCurve(painter, xMap, yMap, polyline);
}


/*!
  \brief Specify an attribute for drawing the curve

  The attributes can be used to modify the drawing style.
  The following attributes are defined:<dl>
  <dt>Fitted</dt>
  <dd>For Lines only. A QwtCurveFitter tries to
      interpolate/smooth the curve, before it is painted.
      Note that curve fitting requires temorary memory
      for calculating coefficients and additional points.
      If painting in Fitted mode is slow it might be better
      to fit the points, before they are passed to PlotCurve.
  </dd>
  <dt>Inverted</dt>
  <dd>For Steps only. Draws a step function
      from the right to the left.</dd></dl>

  \param attribute Curve attribute
  \param on On/Off

  /sa testCurveAttribute(), setCurveFitter()
*/
void PlotCurve::setCurveAttribute(CurveAttribute attribute, bool on)
{
    if ( bool(d_data->attributes & attribute) == on )
        return;

    if ( on )
        d_data->attributes |= attribute;
    else
        d_data->attributes &= ~attribute;

    itemChanged();
}

/*!
    Return the current curve attributes
    \sa setCurveAttribute()
*/
bool PlotCurve::testCurveAttribute(CurveAttribute attribute) const
{
    return d_data->attributes & attribute;
}

/*!
  Assign the curve type

  <dt>PlotCurve::Yfx
  <dd>Draws y as a function of x (the default). The
      baseline is interpreted as a horizontal line
      with y = baseline().</dd>
  <dt>PlotCurve::Xfy
  <dd>Draws x as a function of y. The baseline is
      interpreted as a vertical line with x = baseline().</dd>

  The baseline is used for aligning the sticks, or
  filling the curve with a brush.

  \sa curveType()
*/
void PlotCurve::setCurveType(CurveType curveType)
{
    if ( d_data->curveType != curveType )
    {
        d_data->curveType = curveType;
        itemChanged();
    }
}

/*!
   Return the curve type
   \sa setCurveType()
*/
PlotCurve::CurveType PlotCurve::curveType() const
{
    return d_data->curveType;
}

void PlotCurve::setCurveFitter(QwtCurveFitter *curveFitter)
{
    delete d_data->curveFitter;
    d_data->curveFitter = curveFitter;

    itemChanged();
}

QwtCurveFitter *PlotCurve::curveFitter() const
{
    return d_data->curveFitter;
}

/*!
  Fill the area between the curve and the baseline with
  the curve brush

  \param painter Painter
  \param xMap x map
  \param yMap y map
  \param pa Polygon

  \sa setBrush(), setBaseline(), setCurveType()
*/

void PlotCurve::fillCurve(QPainter *painter,
    const QwtScaleMap &xMap, const QwtScaleMap &yMap,
    QwtPolygon &pa) const
{
    if ( d_data->brush.style() == Qt::NoBrush )
        return;

    closePolyline(xMap, yMap, pa);
    if ( pa.count() <= 2 ) // a line can't be filled
        return;

    QBrush b = d_data->brush;
    if ( !b.color().isValid() )
        b.setColor(d_data->pen.color());

    painter->save();

    painter->setPen(QPen(Qt::NoPen));
    painter->setBrush(b);

    QwtPainter::drawPolygon(painter, pa);

    painter->restore();
}

/*!
  \brief Complete a polygon to be a closed polygon
         including the area between the original polygon
         and the baseline.
  \param xMap X map
  \param yMap Y map
  \param pa Polygon to be completed
*/

void PlotCurve::closePolyline(
    const QwtScaleMap &xMap, const QwtScaleMap &yMap,
    QwtPolygon &pa) const
{
    const int sz = pa.size();
    if ( sz < 2 )
        return;

    pa.resize(sz + 2);

    if ( d_data->curveType == PlotCurve::Xfy )
    {
        pa.setPoint(sz,
            xMap.transform(d_data->reference), pa.point(sz - 1).y());
        pa.setPoint(sz + 1,
            xMap.transform(d_data->reference), pa.point(0).y());
    }
    else
    {
        pa.setPoint(sz,
            pa.point(sz - 1).x(), yMap.transform(d_data->reference));
        pa.setPoint(pa.size() - 1,
            pa.point(0).x(), yMap.transform(d_data->reference));
    }
}

/*!
  \brief Draw symbols
  \param painter Painter
  \param symbol Curve symbol
  \param xMap x map
  \param yMap y map
  \param from index of the first point to be painted
  \param to index of the last point to be painted

  \sa setSymbol(), draw(), drawCurve()
*/
void PlotCurve::drawSymbols(QPainter *painter, const QwtSymbol &symbol,
    const QwtScaleMap &xMap, const QwtScaleMap &yMap,
    int from, int to) const
{
    painter->setBrush(symbol.brush());
    painter->setPen(symbol.pen());

    QRect rect;
    rect.setSize(QwtPainter::metricsMap().screenToLayout(symbol.size()));

    if ( to > from && d_data->paintAttributes & PaintFiltered )
    {
        const QRect window = painter->window();
        if ( window.isEmpty() )
            return;

        PrivateData::PixelMatrix pixelMatrix(window);

        for (int i = from; i <= to; i++)
        {
            const QPoint pi( xMap.transform(x(i)),
                yMap.transform(y(i)) );

            if ( pixelMatrix.testPixel(pi) )
            {
                rect.moveCenter(pi);
                symbol.draw(painter, rect);
            }
        }
    }
    else
    {
        for (int i = from; i <= to; i++)
        {
            const int xi = xMap.transform(x(i));
            const int yi = yMap.transform(y(i));

            rect.moveCenter(QPoint(xi, yi));
            symbol.draw(painter, rect);
        }
    }
}

/*!
  \brief Set the value of the baseline

  The baseline is needed for filling the curve with a brush or
  the Sticks drawing style.
  The default value is 0.0. The interpretation
  of the baseline depends on the CurveType. With PlotCurve::Yfx,
  the baseline is interpreted as a horizontal line at y = baseline(),
  with PlotCurve::Yfy, it is interpreted as a vertical line at
  x = baseline().
  \param reference baseline
  \sa baseline(), setBrush(), setStyle(), setCurveType()
*/
void PlotCurve::setBaseline(double reference)
{
    if ( d_data->reference != reference )
    {
        d_data->reference = reference;
        itemChanged();
    }
}

/*!
    Return the value of the baseline
    \sa setBaseline
*/
double PlotCurve::baseline() const
{
    return d_data->reference;
}

/*!
  Return the size of the data arrays
  \sa setData()
*/
int PlotCurve::dataSize() const
{
    return d_xy->size();
}

int PlotCurve::closestPoint(const QPoint &pos, double *dist) const
{
    if ( plot() == NULL || dataSize() <= 0 )
        return -1;

    const QwtScaleMap xMap = plot()->canvasMap(xAxis());
    const QwtScaleMap yMap = plot()->canvasMap(yAxis());

    int index = -1;
    double dmin = 1.0e10;

    for (int i=0; i < dataSize(); i++)
    {
        const double cx = xMap.xTransform(x(i)) - pos.x();
        const double cy = yMap.xTransform(y(i)) - pos.y();

        const double f = qwtSqr(cx) + qwtSqr(cy);
        if (f < dmin)
        {
            index = i;
            dmin = f;
        }
    }
    if ( dist )
        *dist = sqrt(dmin);

    return index;
}

void PlotCurve::updateLegend(QwtLegend *legend) const
{
    if ( !legend )
        return;

    QwtPlotItem::updateLegend(legend);

    QWidget *widget = legend->find(this);
    if ( !widget || !widget->inherits("QwtLegendItem") )
        return;

    QwtLegendItem *legendItem = (QwtLegendItem *)widget;

#if QT_VERSION < 0x040000
    const bool doUpdate = legendItem->isUpdatesEnabled();
#else
    const bool doUpdate = legendItem->updatesEnabled();
#endif
    legendItem->setUpdatesEnabled(false);

    const int policy = legend->displayPolicy();

    if (policy == QwtLegend::FixedIdentifier)
    {
        int mode = legend->identifierMode();

        if (mode & QwtLegendItem::ShowLine)
            legendItem->setCurvePen(pen());

        if (mode & QwtLegendItem::ShowSymbol)
            legendItem->setSymbol(symbol());

        if (mode & QwtLegendItem::ShowText)
            legendItem->setText(title());
        else
            legendItem->setText(QwtText());

        legendItem->setIdentifierMode(mode);
    }
    else if (policy == QwtLegend::AutoIdentifier)
    {
        int mode = 0;

        if (PlotCurve::NoCurve != style())
        {
            legendItem->setCurvePen(pen());
            mode |= QwtLegendItem::ShowLine;
        }
        if (QwtSymbol::NoSymbol != symbol().style())
        {
            legendItem->setSymbol(symbol());
            mode |= QwtLegendItem::ShowSymbol;
        }
        if ( !title().isEmpty() )
        {
            legendItem->setText(title());
            mode |= QwtLegendItem::ShowText;
        }
        else
        {
            legendItem->setText(QwtText());
        }
        legendItem->setIdentifierMode(mode);
    }

    legendItem->setUpdatesEnabled(doUpdate);
    legendItem->update();
}

} // end namespace rbm

