/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Implementation of the QGraphicsGrid class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "QGraphicsGrid.h"
#include <QPainter>

/**
 * Creates a grid using a QRectF object and a default pen and no ticks.
 * \param rect A rectangle specifying the size of the grid.
 * \param xCells The number of grid cells in x-direction.
 * \param yCells The number of grid cells in y-direction.
 * \param parent The parent QGraphicsItem.
 */
QGraphicsGrid::QGraphicsGrid(QRectF rect, int xCells, int yCells,
                             QGraphicsItem* parent) : QGraphicsItem(parent)
{
    _numberOfXCells = xCells;
    _numberOfYCells = yCells;
    _bounds        = rect;
    _showTicks     = false;

    initDefaultPens();
}

/**
 * Creates a grid by specifying its bounding rectangle using a default pen and no ticks.
 * \param x X-coordinate for the lower left corner of the bounding rectangle for the grid.
 * \param y Y-coordinate for the lower left corner of the bounding rectangle for the grid.
 * \param width Width of the bounding rectangle.
 * \param height Height of the bounding rectangle.
 * \param xCells The number of grid cells in x-direction.
 * \param yCells The number of grid cells in y-direction.
 * \param parent The parent QGraphicsItem.
 */
QGraphicsGrid::QGraphicsGrid(int x,
                             int y,
                             int width,
                             int height,
                             int xCells,
                             int yCells,
                             QGraphicsItem* parent) : QGraphicsItem(parent)
{
    _numberOfXCells = xCells;
    _numberOfYCells = yCells;
    _bounds         = QRectF(x,y,width,height);
    _showTicks      = false;

    initDefaultPens();
}

/**
 * Creates a grid using a QRectF object and a default pen and no ticks.
 * \param rect A rectangle specifying the size of the grid.
 * \param xCells The number of grid cells in x-direction.
 * \param yCells The number of grid cells in y-direction.
 * \param ticks Specifies if ticks are displayed for the grid.
 * \param pen The pen for drawing the grid.
 * \param parent The parent QGraphicsItem.
 */
QGraphicsGrid::QGraphicsGrid(QRectF rect,
                             int xCells,
                             int yCells,
                             bool ticks,
                             QPen pen,
                             QGraphicsItem* parent) : QGraphicsItem(parent)
{
    _numberOfXCells = xCells;
    _numberOfYCells = yCells;
    _bounds         = rect;
    _showTicks      = ticks;

    _outside = pen;
    _outside.setCosmetic(true);

    _inside  = pen;
    QColor iColour = pen.color();
    iColour.setAlpha(125);
    _inside.setColor(iColour);
    _inside.setStyle(Qt::DotLine);
    _inside.setCosmetic(true);
}

/**
 * Creates a grid by specifying its bounding rectangle using a default pen and no ticks.
 * \param x X-coordinate for the lower left corner of the bounding rectangle for the grid.
 * \param y Y-coordinate for the lower left corner of the bounding rectangle for the grid.
 * \param width Width of the bounding rectangle.
 * \param height Height of the bounding rectangle.
 * \param xCells The number of grid cells in x-direction.
 * \param yCells The number of grid cells in y-direction.
 * \param ticks Specifies if ticks are displayed for the grid.
 * \param pen The pen for drawing the grid.
 * \param parent The parent QGraphicsItem.
 */
QGraphicsGrid::QGraphicsGrid(int x,
                             int y,
                             int width,
                             int height,
                             int xCells,
                             int yCells,
                             bool ticks,
                             QPen pen,
                             QGraphicsItem* parent) : QGraphicsItem(parent)
{
    _numberOfXCells = xCells;
    _numberOfYCells = yCells;
    _bounds         = QRectF(x,y,width,height);
    _showTicks      = ticks;

    _outside = pen;
    _outside.setCosmetic(true);

    _inside  = pen;
    QColor iColour = pen.color();
    iColour.setAlpha(125);
    _inside.setColor(iColour);
    _inside.setStyle(Qt::DotLine);
    _inside.setCosmetic(true);
}

QGraphicsGrid::~QGraphicsGrid() = default;

/// Returns the bounding rectangle of the grid.
QRectF QGraphicsGrid::boundingRect() const
{
    return _bounds;
}

/// Defines the default pens.
void QGraphicsGrid::initDefaultPens()
{
    QPen in(Qt::gray, 1, Qt::DotLine, Qt::SquareCap, Qt::RoundJoin);
    QPen out(Qt::black, 1, Qt::SolidLine, Qt::SquareCap, Qt::RoundJoin);
    _inside  = in;
    _outside = out;
    _inside.setCosmetic(true);
    _outside.setCosmetic(true);
}

/// Paints the grid.
void QGraphicsGrid::paint(QPainter* painter,
                          const QStyleOptionGraphicsItem* option,
                          QWidget* widget)
{
    Q_UNUSED (option)
    Q_UNUSED (widget)

    if (!_bounds.isValid())
        return;

    /* draw outside rectangle */
    QBrush brush(Qt::NoBrush);
    painter->setPen(_outside);
    painter->drawRect(_bounds);

    /* draw horizontal lines */
    for (int i = 0; i <= _numberOfXCells; ++i)
    {
        auto x = static_cast<int>(
            _bounds.left() + (i * (_bounds.width() - 1) / _numberOfXCells));

        if (i > 0 && i < _numberOfXCells)
        {
            painter->setPen(_inside);
            painter->drawLine(x, (int)_bounds.top(), x, (int)_bounds.bottom());
        }

        /* draw ticks on x-axis */
        if (_showTicks)
        {
            //double label = bounds.left() + (i * bounds.width() / numberOfXCells);
            painter->setPen(_outside);
            painter->drawLine(x, (int)_bounds.bottom(), x, (int)_bounds.bottom() + 5);
            //painter->drawText(x - margin, bounds.bottom() + 5, 2*margin, 20,
            //                   Qt::AlignHCenter | Qt::AlignTop, QString::number(label));
        }
    }

    /* draw vertical lines */
    for (int j = 0; j <= _numberOfYCells; ++j)
    {
        auto y = static_cast<int>(
            _bounds.bottom() - (j * (_bounds.height() - 1) / _numberOfYCells));

        if (j > 0 && j < _numberOfYCells)
        {
            painter->setPen(_inside);
            painter->drawLine((int)_bounds.left(), y, (int)_bounds.right(), y);
        }

        /* draw ticks on y-axis */
        if (_showTicks)
        {
            //double label = bounds.top() + (j * bounds.height() / numberOfYCells);
            painter->setPen(_outside);
            painter->drawLine((int)_bounds.left() - 5, y, (int)_bounds.left(), y);
            //painter->drawText(bounds.left() - margin, y - 10, margin - 5, 20,
            //                   Qt::AlignRight | Qt::AlignVCenter, QString::number(label));
        }
    }
}

/// Sets the number of cells in x direction.
void QGraphicsGrid::setNumberOfXCells(int xCells)
{
    _numberOfXCells = xCells;
}

/// Sets the number of cells in y direction.
void QGraphicsGrid::setNumberOfYCells(int yCells)
{
    _numberOfYCells = yCells;
}

/// Sets the bounding rectangle of the grid.
void QGraphicsGrid::setRect(int x, int y, int width, int height)
{
    _bounds = QRectF(x,y,width,height);
}
