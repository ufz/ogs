// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
    explicit QNonScalableGraphicsTextItem(QGraphicsItem* parent = nullptr);
    explicit QNonScalableGraphicsTextItem(const QString& text,
                                          QGraphicsItem* parent = nullptr);
    ~QNonScalableGraphicsTextItem() override;

    void paint(QPainter* painter,
               const QStyleOptionGraphicsItem* option,
               QWidget* widget) override;
    QRectF boundingRect() const override;
};
