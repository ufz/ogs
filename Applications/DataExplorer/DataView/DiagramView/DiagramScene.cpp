/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Implementation of the DiagramScene class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DiagramScene.h"
#include <limits>
#include <cmath>

// default size of a new window
const float DEFAULTX = 500.0;
const float DEFAULTY = 300.0;

/**
 * Creates a new scene. Since no data points are given some default
 * values are used for constructing all the necessary objects.
 */
DiagramScene::DiagramScene(QObject* parent) : QGraphicsScene(parent)
{
    bounds_.setRect(0,0,1,1);
    initialize();
}

/**
 * Creates a new scene.
 * \param list includes all necessary information of the graph to display.
 * \param parent The parent QObject.
 */
DiagramScene::DiagramScene(DiagramList* list, QObject* parent) : QGraphicsScene(parent)
{
    setDiagramBoundaries(list);
    initialize();
}

DiagramScene::~DiagramScene()
{
    delete grid_;
    delete xAxis_;
    delete yAxis_;
    delete xLabel_;
    delete yLabel_;
    delete xUnit_;
    delete yUnit_;
    for (auto& graphCaption : graphCaptions_)
    {
        delete graphCaption;
    }
    graphCaptions_.clear();
    for (auto& graph : graphs_)
    {
        delete graph;
    }
    graphs_.clear();
    for (auto& text : xTicksText_)
    {
        delete text;
    }
    xTicksText_.clear();
    for (auto& text : yTicksText_)
    {
        delete text;
    }
    yTicksText_.clear();
    for (auto& list : lists_)
    {
        delete list;
    }
    lists_.clear();
}

/// Adds an arrow object to the diagram which might be used as a coordinate axis, etc.
QArrow* DiagramScene::addArrow(qreal length, qreal angle, QPen& pen)
{
    auto* arrow = new QArrow(length, angle, 8, 5, pen);
    addItem(arrow);
    return arrow;
}

/// Adds a caption for a graph beneath the actual diagram.
void DiagramScene::addCaption(const QString &name, QPen &pen)
{
    auto* caption = new QGraphicsItemGroup(nullptr);
    QGraphicsLineItem* l = addLine(0,0,100,0,pen);
    QGraphicsTextItem* t = addText(name);
    l->setPos(0,0);
    t->setPos(110, -(t->boundingRect()).height() / 2);
    caption->addToGroup(l);
    caption->addToGroup(t);
    caption->setFlag(QGraphicsItem::ItemIgnoresTransformations, true);

    graphCaptions_.push_back(caption);
    addItem(graphCaptions_[graphCaptions_.size() - 1]);
}

/// Adds a graph to the scene, including all data points and meta-information.
void DiagramScene::addGraph(DiagramList* list)
{
    setDiagramBoundaries(list);
    adjustScaling();
    xLabel_->setPlainText(list->getXLabel());
    yLabel_->setPlainText(list->getYLabel());
    xUnit_->setPlainText(list->getXUnit());
    yUnit_->setPlainText(list->getYUnit());

    clearGrid();
    constructGrid();

    lists_.push_back(list);
    for (auto& list : lists_)
    {
        drawGraph(list);
    }

    update();
}

/// Adds a grid-object to the scene
QGraphicsGrid* DiagramScene::addGrid(const QRectF &rect, int xTicks, int yTicks, const QPen &pen)
{
    QGraphicsGrid* g = new QGraphicsGrid(rect, xTicks, yTicks, true, pen);
    addItem(g);
    return g;
}

/// Adds a non-scalable text object to the scene
QNonScalableGraphicsTextItem* DiagramScene::addNonScalableText(const QString &text,
                                                               const QFont &font)
{
    auto* item = new QNonScalableGraphicsTextItem(text);
    item->setFont(font);
    addItem(item);
    return item;
}

/// Resizes a given axis to "nice" dimensions and calculates an adequate number of ticks to be placed on it
void DiagramScene::adjustAxis(qreal& min, qreal& max, int& numberOfTicks)
{
    const int MinTicks = 4;
    double grossStep = (max - min) / MinTicks;
    double step = pow(10.0, std::floor(log10(grossStep)));
    if (5 * step < grossStep)
    {
        step *= 5;
    }
    else if (2 * step < grossStep)
    {
        step *= 2;
    }
    numberOfTicks = int(ceil(max / step) - std::floor(min / step));
    if (numberOfTicks < MinTicks)
    {
        numberOfTicks = MinTicks;
    }
    min = std::floor(min / step) * step;
    max = ceil(max / step) * step;
}

///Calculates scaling factors to set coordinate system and graphs to default window size
void DiagramScene::adjustScaling()
{
    if ( (unscaledBounds_.width() > 0) && (unscaledBounds_.height() > 0))
    {
        scaleX_ = DEFAULTX / static_cast<float>(unscaledBounds_.width());
        scaleY_ = DEFAULTY / static_cast<float>(unscaledBounds_.height());
    }
}

/// Destroys the grid object (coordinate system) when a new graph is added.
void DiagramScene::clearGrid()
{
    if (!lists_.isEmpty())
    {
        removeItem(grid_);

        for (auto& text : xTicksText_)
        {
            removeItem(text);
        }
        for (auto& text : yTicksText_)
        {
            removeItem(text);
        }
        for (auto& graph : graphs_)
        {
            removeItem(graph);
        }
        for (auto& graphCaption : graphCaptions_)
        {
            removeItem(graphCaption);
        }

        xTicksText_.clear();
        yTicksText_.clear();
        graphs_.clear();
        graphCaptions_.clear();
    }
}

/// Adjusts the underlying grid based on the graphs that are displayed in the diagram
void DiagramScene::constructGrid()
{
    // be very careful with scaling parameters here!
    int numXTicks;
    int numYTicks;
    qreal xMin = unscaledBounds_.left();
    qreal yMin = unscaledBounds_.top();
    qreal xMax = unscaledBounds_.right();
    qreal yMax = unscaledBounds_.bottom();

    adjustAxis(xMin, xMax, numXTicks);
    adjustAxis(yMin, yMax, numYTicks);

    // adjust boundaries of coordinate system according to scaling
    bounds_.setRect(    xMin * scaleX_,
                        yMin * scaleY_,
                        (xMax - xMin) * scaleX_,
                        (yMax - yMin) * scaleY_
                        );

    QPen pen(Qt::black, 1, Qt::SolidLine, Qt::SquareCap, Qt::RoundJoin);
    grid_ = addGrid(bounds_, numXTicks, numYTicks, pen);

    if (startDate_ == QDateTime())
    {
        for (int i = 0; i <= numXTicks; ++i)
        {
            auto x =
                static_cast<int>(bounds_.left() / scaleX_ +
                                 (i * (bounds_.width() / scaleX_) / numXTicks));
            xTicksText_.push_back(addNonScalableText(QString::number(x)));
            xTicksText_.last()->setPos(x * scaleX_, bounds_.bottom() + 15);
        }
    }
    else
    {
        for (int i = 0; i <= numXTicks; ++i)
        {
            auto x =
                static_cast<int>(bounds_.left() / scaleX_ +
                                 (i * (bounds_.width() / scaleX_) / numXTicks));
            QDateTime currentDate = startDate_.addSecs(x);
            xTicksText_.push_back(
                addNonScalableText(currentDate.toString("dd.MM.yyyy")));
            xTicksText_.last()->setPos(x * scaleX_, bounds_.bottom() + 15);
        }
    }

    for (int j = 0; j <= numYTicks; ++j)
    {
        qreal y = bounds_.bottom() / scaleY_ -
                  (j * (bounds_.height() / scaleY_) / numYTicks);
        qreal label = bounds_.top() / scaleY_ +
                      (j * (bounds_.height() / scaleY_) / numYTicks);
        yTicksText_.push_back(addNonScalableText(QString::number(label)));
        yTicksText_.last()->setPos(bounds_.left() - MARGIN / 2, y * scaleY_);
    }
}

/// Plots the graph.
void DiagramScene::drawGraph(DiagramList* list)
{
    QPainterPath path;

    if (list->getPath(path, scaleX_, scaleY_))
    {
        QPen pen(list->getColor(), 2, Qt::SolidLine, Qt::SquareCap, Qt::RoundJoin);
        pen.setCosmetic(true);
        graphs_.push_back(addPath(path, pen));
        addCaption(list->getName(), pen);

        int last = graphs_.size() - 1;

        /**
         * For correct display the graph needs to be flipped vertically and then
         * translated back to its original position
         */
        int verticalShift =
                static_cast<int>(2 *
                                 (list->minYValue() *
                                  scaleY_) + (graphs_[last]->boundingRect()).height());
        graphs_[last]->setTransform(QTransform(QMatrix(1,0,0,-1,0,verticalShift)));
    }
}

/// Returns the y-value at which the x-axis should cross the y-axis.
/// This value is zero if minYValue<0<maxYValue and minYValue otherwise.
int DiagramScene::getXAxisOffset()
{
    return (bounds_.top() <= 0 && bounds_.bottom() > 0)
               ? static_cast<int>(bounds_.bottom() + bounds_.top())
               : static_cast<int>(bounds_.bottom());
}

/// Returns the x-value at which the y-axis should cross the x-axis.
/// This value is zero if minXValue<0<maxXValue and minXValue otherwise.
int DiagramScene::getYAxisOffset()
{
    return (bounds_.left() <= 0 && bounds_.right() > 0)
               ? 0
               : static_cast<int>(bounds_.left());
}

/// Initialises the coordinate axes, adds labels and/or units to the axes,
/// calculates axes-lengths, offsets, etc.
void DiagramScene::initialize()
{
    QPen pen(Qt::black, 1, Qt::SolidLine, Qt::SquareCap, Qt::RoundJoin);
    pen.setCosmetic(true);

    setXAxis(addArrow(bounds_.width(),  0, pen));
    setYAxis(addArrow(bounds_.height(), -90, pen));
    xLabel_ = addNonScalableText(" ");
    yLabel_ = addNonScalableText(" ");
    yLabel_->setRotation(-90);

    xUnit_ = addNonScalableText(" ");
    yUnit_ = addNonScalableText(" ");

    update();
}

/// Updates the (unscaled) boundaries of the visible coordinate system when a new
/// list is added (boundaries are rescaled in the constructGrid-method
void DiagramScene::setDiagramBoundaries(DiagramList* list)
{
    if (!lists_.isEmpty())
    {
        if (list->minXValue() < unscaledBounds_.left())
        {
            unscaledBounds_.setLeft(list->minXValue());
        }
        if (list->minYValue() < unscaledBounds_.top())
        {
            unscaledBounds_.setTop(list->minYValue());
        }
        if (list->maxXValue() > unscaledBounds_.right())
        {
            unscaledBounds_.setRight(list->maxXValue());
        }
        if (list->maxYValue() > unscaledBounds_.bottom())
        {
            unscaledBounds_.setBottom(list->maxYValue());
        }
        if (startDate_ > list->getStartDate())
        {
            startDate_ = list->getStartDate();
        }
    }
    else
    {
        unscaledBounds_.setRect(list->minXValue(), list->minYValue(),
                                list->maxXValue() - list->minXValue(),
                                list->maxYValue() - list->minYValue());
        startDate_ = list->getStartDate();
    }
}

/**
 * Updates the scene at the start and everytime new data points
 * are added. Specifically, objects on the scene are assigned
 * their position in the new coordinate system and are resized
 * if necessary.
 */
void DiagramScene::update()
{
    xAxis_->setPos(bounds_.left(),getXAxisOffset());
    yAxis_->setPos(getYAxisOffset(),bounds_.bottom());
    xAxis_->setLength(bounds_.width());
    yAxis_->setLength(bounds_.height());

    xLabel_->setPos( bounds_.left() + bounds_.width() / 2, bounds_.bottom() + 1.5 * MARGIN );
    yLabel_->setPos( bounds_.left() - 1.5 * MARGIN, bounds_.top() + bounds_.height() / 2 );

    xUnit_->setPos( bounds_.right(), bounds_.bottom() + 1.2 * MARGIN);
    yUnit_->setPos( bounds_.left(), bounds_.top() - 0.5 * MARGIN);

    /* update graphs and their captions */
    QRectF rect;
    for (int i = 0; i < graphs_.size(); i++)
    {
        rect = graphs_[i]->boundingRect();
        auto offset = static_cast<int>(fabs(rect.bottom() - bounds_.bottom()) -
                                       fabs(rect.top() - bounds_.top()));
        graphs_[i]->setPos(0, offset);

        rect = itemsBoundingRect();
        graphCaptions_[i]->setPos(bounds_.left(),rect.bottom() + 10);
    }
}
