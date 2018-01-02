/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Definition of the QGraphicsGrid class.
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

/**
 * \brief A 2D carthesian grid as a QGraphicsItem.
 *
 * A 2D carthesian grid as a QGraphicsItem. The size of the grid cells is constant but can be anisotroph.
 */
class QGraphicsGrid : public QGraphicsItem
{
public:
    QGraphicsGrid(QRectF rect,
                  int xCells,
                  int yCells,
                  QGraphicsItem* parent = nullptr);
    QGraphicsGrid(int x,
                  int y,
                  int width,
                  int height,
                  int xCells,
                  int yCells,
                  QGraphicsItem* parent = nullptr);
    QGraphicsGrid(QRectF rect,
                  int xCells,
                  int yCells,
                  bool ticks,
                  QPen pen,
                  QGraphicsItem* parent = nullptr);
    QGraphicsGrid(int x,
                  int y,
                  int width,
                  int height,
                  int xCells,
                  int yCells,
                  bool ticks,
                  QPen pen,
                  QGraphicsItem* parent = nullptr);
    ~QGraphicsGrid() override;

    QRectF boundingRect() const override;
    void paint(QPainter* painter,
               const QStyleOptionGraphicsItem* option,
               QWidget* widget) override;
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
