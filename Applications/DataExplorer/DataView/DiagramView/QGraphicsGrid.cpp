/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Implementation of the QGraphicsGrid class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
    numberOfXCells_ = xCells;
    numberOfYCells_ = yCells;
    bounds_        = rect;
    showTicks_     = false;

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
    numberOfXCells_ = xCells;
    numberOfYCells_ = yCells;
    bounds_         = QRectF(x,y,width,height);
    showTicks_      = false;

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
    numberOfXCells_ = xCells;
    numberOfYCells_ = yCells;
    bounds_         = rect;
    showTicks_      = ticks;

    outside_ = pen;
    outside_.setCosmetic(true);

    inside_  = pen;
    QColor iColour = pen.color();
    iColour.setAlpha(125);
    inside_.setColor(iColour);
    inside_.setStyle(Qt::DotLine);
    inside_.setCosmetic(true);
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
    numberOfXCells_ = xCells;
    numberOfYCells_ = yCells;
    bounds_         = QRectF(x,y,width,height);
    showTicks_      = ticks;

    outside_ = pen;
    outside_.setCosmetic(true);

    inside_  = pen;
    QColor iColour = pen.color();
    iColour.setAlpha(125);
    inside_.setColor(iColour);
    inside_.setStyle(Qt::DotLine);
    inside_.setCosmetic(true);
}

QGraphicsGrid::~QGraphicsGrid() = default;

/// Returns the bounding rectangle of the grid.
QRectF QGraphicsGrid::boundingRect() const
{
    return bounds_;
}

/// Defines the default pens.
void QGraphicsGrid::initDefaultPens()
{
    QPen in(Qt::gray, 1, Qt::DotLine, Qt::SquareCap, Qt::RoundJoin);
    QPen out(Qt::black, 1, Qt::SolidLine, Qt::SquareCap, Qt::RoundJoin);
    inside_  = in;
    outside_ = out;
    inside_.setCosmetic(true);
    outside_.setCosmetic(true);
}

/// Paints the grid.
void QGraphicsGrid::paint(QPainter* painter,
                          const QStyleOptionGraphicsItem* option,
                          QWidget* widget)
{
    Q_UNUSED (option)
    Q_UNUSED (widget)

    if (!bounds_.isValid())
    {
        return;
    }

    /* draw outside rectangle */
    QBrush brush(Qt::NoBrush);
    painter->setPen(outside_);
    painter->drawRect(bounds_);

    /* draw horizontal lines */
    for (int i = 0; i <= numberOfXCells_; ++i)
    {
        auto x = static_cast<int>(
            bounds_.left() + (i * (bounds_.width() - 1) / numberOfXCells_));

        if (i > 0 && i < numberOfXCells_)
        {
            painter->setPen(inside_);
            painter->drawLine(x, static_cast<int>(bounds_.top()), x,
                              static_cast<int>(bounds_.bottom()));
        }

        /* draw ticks on x-axis */
        if (showTicks_)
        {
            //double label = bounds.left() + (i * bounds.width() / numberOfXCells);
            painter->setPen(outside_);
            painter->drawLine(x, static_cast<int>(bounds_.bottom()), x,
                              static_cast<int>(bounds_.bottom()) + 5);
            //painter->drawText(x - margin, bounds.bottom() + 5, 2*margin, 20,
            //                   Qt::AlignHCenter | Qt::AlignTop, QString::number(label));
        }
    }

    /* draw vertical lines */
    for (int j = 0; j <= numberOfYCells_; ++j)
    {
        auto y = static_cast<int>(
            bounds_.bottom() - (j * (bounds_.height() - 1) / numberOfYCells_));

        if (j > 0 && j < numberOfYCells_)
        {
            painter->setPen(inside_);
            painter->drawLine(static_cast<int>(bounds_.left()), y,
                              static_cast<int>(bounds_.right()), y);
        }

        /* draw ticks on y-axis */
        if (showTicks_)
        {
            //double label = bounds.top() + (j * bounds.height() / numberOfYCells);
            painter->setPen(outside_);
            painter->drawLine(static_cast<int>(bounds_.left()) - 5, y,
                              static_cast<int>(bounds_.left()), y);
            //painter->drawText(bounds.left() - margin, y - 10, margin - 5, 20,
            //                   Qt::AlignRight | Qt::AlignVCenter, QString::number(label));
        }
    }
}

/// Sets the number of cells in x direction.
void QGraphicsGrid::setNumberOfXCells(int xCells)
{
    numberOfXCells_ = xCells;
}

/// Sets the number of cells in y direction.
void QGraphicsGrid::setNumberOfYCells(int yCells)
{
    numberOfYCells_ = yCells;
}

/// Sets the bounding rectangle of the grid.
void QGraphicsGrid::setRect(int x, int y, int width, int height)
{
    bounds_ = QRectF(x,y,width,height);
}
