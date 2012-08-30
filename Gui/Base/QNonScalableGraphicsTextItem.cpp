/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file QNonScalableGraphicsTextItem.cpp
 *
 * Created on by Karsten Rink
 */

#include "QNonScalableGraphicsTextItem.h"
#include <QPainter>

/// Constructor using a QGraphicsTextItem.
QNonScalableGraphicsTextItem::QNonScalableGraphicsTextItem(QGraphicsItem* parent) :
	QGraphicsTextItem(parent)
{
	setAcceptDrops(true);
	setAcceptHoverEvents(true);
	setFlag(QGraphicsItem::ItemIgnoresTransformations, true);
}

/// Constructor using a QString.
QNonScalableGraphicsTextItem::QNonScalableGraphicsTextItem(const QString & text,
                                                           QGraphicsItem* parent) :
	QGraphicsTextItem(parent)
{
	if (!text.isEmpty())
		setPlainText(text);
	setAcceptDrops(true);
	setAcceptHoverEvents(true);
	setFlag(QGraphicsItem::ItemIgnoresTransformations, true);
}

QNonScalableGraphicsTextItem::~QNonScalableGraphicsTextItem()
{
}

/// Paints the text item.
void QNonScalableGraphicsTextItem::paint(QPainter* painter,
                                         const QStyleOptionGraphicsItem* option,
                                         QWidget* widget)
{
	//painter->drawRect(boundingRect());
	QRectF rect = boundingRect();
	painter->translate(-rect.width() / 2, -rect.height() / 2);
	QGraphicsTextItem::paint(painter, option, widget);
}

/// Returns the bounding rectangle of the text item.
QRectF QNonScalableGraphicsTextItem::boundingRect() const
{
	QRectF rect = QGraphicsTextItem::boundingRect();
	return rect; //QRectF(rect.x()-rect.width()/2, rect.y()-rect.height()/2,rect.width(), rect.height());
}
