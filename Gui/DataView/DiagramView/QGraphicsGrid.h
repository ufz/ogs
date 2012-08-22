/**
 * \file QGraphicsGrid.h
 * KR Initial implementation
 */

#ifndef QGRAPHICSGRID_H
#define QGRAPHICSGRID_H

#include <QGraphicsItem>
#include <QPen>

/**
 * \brief A 2D carthesian grid as a QGraphicsItem.
 *
 * A 2D carthesian grid as a QGraphicsItem. The size of the grid cells is constant but can be anisotroph.
 */
class QGraphicsGrid : public QGraphicsItem
{
public:
	QGraphicsGrid(QRectF rect, int xCells, int yCells, QGraphicsItem* parent = 0);
	QGraphicsGrid(int x,
	              int y,
	              int width,
	              int height,
	              int xCells,
	              int yCells,
	              QGraphicsItem* parent = 0);
	QGraphicsGrid(QRectF rect,
	              int xCells,
	              int yCells,
	              bool ticks,
	              QPen pen,
	              QGraphicsItem* parent = 0);
	QGraphicsGrid(int x,
	              int y,
	              int width,
	              int height,
	              int xCells,
	              int yCells,
	              bool ticks,
	              QPen pen,
	              QGraphicsItem* parent = 0);
	~QGraphicsGrid();

	QRectF boundingRect() const;
	void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget);
	void setNumberOfXCells(int xCells);
	void setNumberOfYCells(int yCells);
	void setRect(QRectF rect);
	void setRect(int x, int y, int width, int height);

private:
	void initDefaultPens();

	QPen _inside;
	QPen _outside;
	QRectF _bounds;
	int _numberOfXCells;
	int _numberOfYCells;
	bool _showTicks;
};

#endif //QGRAPHICSGRID_H
