/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Definition of the QArrow class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef QARROW_H
#define QARROW_H

#include <QGraphicsItem>
#include <QPen>

const double PI = 3.14159265;

/**
 * \brief An arrow as a QGraphicsObject
 */
class QArrow : public QGraphicsItem
{
public:
    QArrow(float l, float a, float hl, float hw, QPen &pen, QGraphicsItem* parent = 0);
    QArrow(float l, float a, QPen &pen, QGraphicsItem* parent = 0);
    ~QArrow();

    double getLength();
    double getAngle();
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget);
    QRectF boundingRect() const;
    void setAngle(double a);
    void setLength(double l);

private:
    double calcCos(double angle);
    double calcSin(double angle);

    float _arrowLength;
    float _arrowAngle;
    float _headLength;
    float _headWidth;
    QPen _arrowPen;
};

#endif //QARROW_H
