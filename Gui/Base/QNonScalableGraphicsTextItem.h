/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file QNonScalableGraphicsTextItem.h
 *
 * Created on by Karsten Rink
 */

#ifndef QNONSCALABLETEXTITEM_H
#define QNONSCALABLETEXTITEM_H

#include <QGraphicsTextItem>

/**
 * \brief A QGraphicsTextItem that will ignore all geometric transformations.
 *
 * A QGraphicsTextItem that will ignore all geometric transformations to the underlying QGraphicsView/QGraphicsScene (in particular, it will not be scaled).
 */
class QNonScalableGraphicsTextItem : public QGraphicsTextItem
{
public:
	QNonScalableGraphicsTextItem(QGraphicsItem* parent = 0);
	QNonScalableGraphicsTextItem(const QString &text, QGraphicsItem* parent = 0);
	~QNonScalableGraphicsTextItem();

	void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget);
	virtual QRectF boundingRect() const;
};

#endif //QNONSCALABLETEXTITEM_H
