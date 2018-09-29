/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Definition of the QNonScalableGraphicsTextItem class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <QGraphicsTextItem>

/**
 * \brief A QGraphicsTextItem that will ignore all geometric transformations.
 *
 * A QGraphicsTextItem that will ignore all geometric transformations to the underlying QGraphicsView/QGraphicsScene (in particular, it will not be scaled).
 */
class QNonScalableGraphicsTextItem : public QGraphicsTextItem
{
public:
    QNonScalableGraphicsTextItem(QGraphicsItem* parent = nullptr);
    QNonScalableGraphicsTextItem(const QString& text,
                                 QGraphicsItem* parent = nullptr);
    ~QNonScalableGraphicsTextItem() override;

    void paint(QPainter* painter,
               const QStyleOptionGraphicsItem* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
};
