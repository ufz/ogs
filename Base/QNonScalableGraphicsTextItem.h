/**
 * \file QNonScalableGraphicsTextItem.h
 * KR Initial implementation
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
	QNonScalableGraphicsTextItem(const QString &text, QGraphicsItem * parent = 0);
	~QNonScalableGraphicsTextItem();

	void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
	QRectF boundingRect();
};

#endif //QNONSCALABLETEXTITEM_H
