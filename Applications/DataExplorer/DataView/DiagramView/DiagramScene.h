/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Definition of the DiagramScene class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef DIAGRAMSCENE_H
#define DIAGRAMSCENE_H

#include "DiagramList.h"
#include "QArrow.h"
#include "QGraphicsGrid.h"
#include "QNonScalableGraphicsTextItem.h"
#include <QGraphicsScene>

class QDateTime;

/**
 * \brief A scene graph for a 2D Diagram including coordinate axes with labels and ticks for one or more plotted graphs.
 */
class DiagramScene : public QGraphicsScene
{
public:
    DiagramScene(QObject* parent = 0);
    DiagramScene(DiagramList* list, QObject* parent = 0);
    ~DiagramScene();

    QArrow* addArrow(float length, float angle, QPen &pen);
    void addGraph(DiagramList* list);
    QGraphicsGrid* addGrid(const QRectF &rect, int xTicks, int yTicks, const QPen &pen);

    static const int MARGIN = 30; /// The margin between the boundary of the scene and the bounding box of all items within the scene

private:
    void addCaption(const QString &name, QPen &pen);
    QNonScalableGraphicsTextItem* addNonScalableText(const QString &text,
                                                     const QFont &font = QFont());
    void adjustAxis(float &min, float &max, int &numberOfTicks);
    void adjustScaling();
    void clearGrid();
    void constructGrid();
    void drawGraph(DiagramList* list);
    int getXAxisOffset();
    int getYAxisOffset();
    void initialize();
    void setDiagramBoundaries(DiagramList* list);

    /// Sets an arrow as x-axis
    void setXAxis(QArrow* arrow) { _xAxis = arrow; }

    /// Sets an arrow as y-axis
    void setYAxis(QArrow* arrow) { _yAxis = arrow; }

    void update();

    QRectF _bounds;
    QRectF _unscaledBounds;
    QVector<DiagramList*> _lists;
    QVector<QGraphicsItemGroup*> _graphCaptions;
    QVector<QGraphicsPathItem*> _graphs;
    QGraphicsGrid* _grid;
    QDateTime _startDate;
    float _scaleX;
    float _scaleY;
    QArrow* _xAxis;
    QArrow* _yAxis;
    QNonScalableGraphicsTextItem* _xLabel;
    QNonScalableGraphicsTextItem* _yLabel;
    QNonScalableGraphicsTextItem* _xUnit;
    QNonScalableGraphicsTextItem* _yUnit;
    QVector<QNonScalableGraphicsTextItem*> _xTicksText;
    QVector<QNonScalableGraphicsTextItem*> _yTicksText;
};

#endif //DIAGRAMSCENE_H
