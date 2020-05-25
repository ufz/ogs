/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Implementation of the QArrow class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "QArrow.h"
#include <QPainter>
#include <math.h>

/**
 * Creates an arrow as a QGraphicItem.
 * \param l Length of the arrow
 * \param a Orientation of the arrow
 * \param hl Length of the arrow head
 * \param hw Width of the arrow head
 * \param pen The pen for drawing the arrow
 * \param parent The parent QGraphicsItem.
 */
QArrow::QArrow(qreal l, qreal a, qreal hl, qreal hw, QPen& pen,
               QGraphicsItem* parent)
    : QGraphicsItem(parent)
{
    arrowLength_ = l;
    arrowAngle_  = a;
    headLength_  = hl;
    headWidth_   = hw;
    arrowPen_    = pen;
}

/**
 * Creates an arrow as a QGraphicItem. Length and width of the arrow head are given by default values.
 * \param l Length of the arrow
 * \param a Orientation of the arrow
 * \param pen The pen for drawing the arrow
 * \param parent The parent QGraphicsItem.
 */
QArrow::QArrow(qreal l, qreal a, QPen& pen, QGraphicsItem* parent)
    : QGraphicsItem(parent)
{
    arrowLength_ = l;
    arrowAngle_  = a;
    headLength_  = 8;   // default headlength
    headWidth_   = 5;   // default headwidth
    arrowPen_    = pen;
}

QArrow::~QArrow() = default;

double QArrow::calcCos(double angle)
{
    return cos (angle * PI / 180);
}

double QArrow::calcSin(double angle)
{
    return sin (angle * PI / 180);
}

/// The bounding box of the arrow
QRectF QArrow::boundingRect() const
{
    double deltaX = cos(arrowAngle_) * arrowLength_;
    double deltaY = sin(arrowAngle_) * arrowLength_;

    return QRectF(0, 0, deltaX, deltaY);
}

/// Returns the length of the arrow.
double QArrow::getLength()
{
    return arrowLength_;
}

/// Returns the orientation of the arrow
double QArrow::getAngle()
{
    return arrowAngle_;
}

/**
 * Overloaded paint-method from QGraphicsItem.
 * Basically it draws a line with an arrowhead consisting of two short lines at the end
 */
void QArrow::paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget)
{
    Q_UNUSED (option)
    Q_UNUSED (widget)

    double ddeltaX    = calcCos(arrowAngle_) * arrowLength_;
    double ddeltaY    = calcSin(arrowAngle_) * arrowLength_;
    double theta     = atan(ddeltaY / ddeltaX);
    double theta2    = (ddeltaX < 0.0) ? (theta + PI) : theta;
    int lengthdeltaX = -static_cast<int>(cos(theta2) * headLength_);
    int lengthdeltaY = -static_cast<int>(sin(theta2) * headLength_);
    auto widthdeltaX = static_cast<int>(sin(theta2) * headWidth_);
    auto widthdeltaY = static_cast<int>(cos(theta2) * headWidth_);
    auto deltaX = static_cast<int>(ddeltaX);
    auto deltaY = static_cast<int>(ddeltaY);
    painter->setPen(arrowPen_);
    painter->drawLine(0, 0, deltaX, deltaY);
    painter->drawLine(deltaX,
                      deltaY,
                      deltaX + lengthdeltaX + widthdeltaX,
                      deltaY + lengthdeltaY - widthdeltaY);
    painter->drawLine(deltaX,
                      deltaY,
                      deltaX + lengthdeltaX - widthdeltaX,
                      deltaY + lengthdeltaY + widthdeltaY);
}

/// Changes orientation of the arrow.
void QArrow::setAngle(qreal a)
{
    arrowAngle_ = a;
}

/// Changes the length of the arrow.
void QArrow::setLength(qreal l)
{
    arrowLength_ = l;
}
