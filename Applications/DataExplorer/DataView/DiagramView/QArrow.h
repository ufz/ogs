/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Definition of the QArrow class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <QGraphicsItem>
#include <QPen>

const double PI = 3.14159265;

/**
 * \brief An arrow as a QGraphicsObject
 */
class QArrow : public QGraphicsItem
{
public:
    QArrow(qreal l, qreal a, qreal hl, qreal hw, QPen& pen,
           QGraphicsItem* parent = nullptr);
    QArrow(qreal l, qreal a, QPen& pen, QGraphicsItem* parent = nullptr);
    ~QArrow() override;

    double getLength();
    double getAngle();
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
    void setAngle(float a);
    void setLength(qreal l);

private:
    double calcCos(double angle);
    double calcSin(double angle);

    qreal _arrowLength;
    qreal _arrowAngle;
    qreal _headLength;
    qreal _headWidth;
    QPen _arrowPen;
};
