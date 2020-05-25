/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Definition of the DiagramScene class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

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
    explicit DiagramScene(QObject* parent = nullptr);
    explicit DiagramScene(DiagramList* list, QObject* parent = nullptr);
    ~DiagramScene() override;

    QArrow* addArrow(qreal length, qreal angle, QPen& pen);
    void addGraph(DiagramList* list);
    QGraphicsGrid* addGrid(const QRectF &rect, int xTicks, int yTicks, const QPen &pen);

    static const int MARGIN = 30; /// The margin between the boundary of the scene and the bounding box of all items within the scene

private:
    void addCaption(const QString &name, QPen &pen);
    QNonScalableGraphicsTextItem* addNonScalableText(const QString &text,
                                                     const QFont &font = QFont());
    void adjustAxis(qreal& min, qreal& max, int& numberOfTicks);
    void adjustScaling();
    void clearGrid();
    void constructGrid();
    void drawGraph(DiagramList* list);
    int getXAxisOffset();
    int getYAxisOffset();
    void initialize();
    void setDiagramBoundaries(DiagramList* list);

    /// Sets an arrow as x-axis
    void setXAxis(QArrow* arrow) { xAxis_ = arrow; }

    /// Sets an arrow as y-axis
    void setYAxis(QArrow* arrow) { yAxis_ = arrow; }

    void update();

    QRectF bounds_;
    QRectF unscaledBounds_;
    QVector<DiagramList*> lists_;
    QVector<QGraphicsItemGroup*> graphCaptions_;
    QVector<QGraphicsPathItem*> graphs_;
    QGraphicsGrid* grid_;
    QDateTime startDate_;
    float scaleX_;
    float scaleY_;
    QArrow* xAxis_;
    QArrow* yAxis_;
    QNonScalableGraphicsTextItem* xLabel_;
    QNonScalableGraphicsTextItem* yLabel_;
    QNonScalableGraphicsTextItem* xUnit_;
    QNonScalableGraphicsTextItem* yUnit_;
    QVector<QNonScalableGraphicsTextItem*> xTicksText_;
    QVector<QNonScalableGraphicsTextItem*> yTicksText_;
};
