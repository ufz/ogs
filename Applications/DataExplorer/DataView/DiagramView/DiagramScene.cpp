/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Implementation of the DiagramScene class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
    _bounds.setRect(0,0,1,1);
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
    delete _grid;
    delete _xAxis;
    delete _yAxis;
    delete _xLabel;
    delete _yLabel;
    delete _xUnit;
    delete _yUnit;
    for (auto& graphCaption : _graphCaptions)
        delete graphCaption;
    _graphCaptions.clear();
    for (auto& graph : _graphs)
        delete graph;
    _graphs.clear();
    for (auto& text : _xTicksText)
        delete text;
    _xTicksText.clear();
    for (auto& text : _yTicksText)
        delete text;
    _yTicksText.clear();
    for (auto& list : _lists)
        delete list;
    _lists.clear();
}

/// Adds an arrow object to the diagram which might be used as a coordinate axis, etc.
QArrow* DiagramScene::addArrow(float length, float angle, QPen &pen)
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

    _graphCaptions.push_back(caption);
    addItem(_graphCaptions[_graphCaptions.size() - 1]);
}

/// Adds a graph to the scene, including all data points and meta-information.
void DiagramScene::addGraph(DiagramList* list)
{
    setDiagramBoundaries(list);
    adjustScaling();
    _xLabel->setPlainText(list->getXLabel());
    _yLabel->setPlainText(list->getYLabel());
    _xUnit->setPlainText(list->getXUnit());
    _yUnit->setPlainText(list->getYUnit());

    clearGrid();
    constructGrid();

    _lists.push_back(list);
    for (auto& list : _lists)
        drawGraph(list);

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
void DiagramScene::adjustAxis(float &min, float &max, int &numberOfTicks)
{
    const int MinTicks = 4;
    double grossStep = (max - min) / MinTicks;
    double step = pow(10.0, std::floor(log10(grossStep)));
    if (5 * step < grossStep)
        step *= 5;
    else if (2 * step < grossStep)
        step *= 2;
    numberOfTicks = int(ceil(max / step) - std::floor(min / step));
    if (numberOfTicks < MinTicks)
        numberOfTicks = MinTicks;
    min = std::floor(min / step) * step;
    max = ceil(max / step) * step;
}

///Calculates scaling factors to set coordinate system and graphs to default window size
void DiagramScene::adjustScaling()
{
    if ( (_unscaledBounds.width() > 0) && (_unscaledBounds.height() > 0))
    {
        _scaleX = DEFAULTX / (_unscaledBounds.width());
        _scaleY = DEFAULTY / (_unscaledBounds.height());
    }
}

/// Destroys the grid object (coordinate system) when a new graph is added.
void DiagramScene::clearGrid()
{
    if (!_lists.isEmpty())
    {
        removeItem(_grid);

        for (auto& text : _xTicksText)
            removeItem(text);
        for (auto& text : _yTicksText)
            removeItem(text);
        for (auto& graph : _graphs)
            removeItem(graph);
        for (auto& graphCaption : _graphCaptions)
            removeItem(graphCaption);

        _xTicksText.clear();
        _yTicksText.clear();
        _graphs.clear();
        _graphCaptions.clear();
    }
}

/// Adjusts the underlying grid based on the graphs that are displayed in the diagram
void DiagramScene::constructGrid()
{
    // be very careful with scaling parameters here!
    int numXTicks, numYTicks;
    float xMin = _unscaledBounds.left();
    float yMin = _unscaledBounds.top();
    float xMax = _unscaledBounds.right();
    float yMax = _unscaledBounds.bottom();

    adjustAxis(xMin, xMax, numXTicks);
    adjustAxis(yMin, yMax, numYTicks);

    // adjust boundaries of coordinate system according to scaling
    _bounds.setRect(    xMin * _scaleX,
                        yMin * _scaleY,
                        (xMax - xMin) * _scaleX,
                        (yMax - yMin) * _scaleY
                        );

    QPen pen(Qt::black, 1, Qt::SolidLine, Qt::SquareCap, Qt::RoundJoin);
    _grid = addGrid(_bounds, numXTicks, numYTicks, pen);

    if (_startDate == QDateTime())
    {
        for (int i = 0; i <= numXTicks; ++i)
        {
            auto x =
                static_cast<int>(_bounds.left() / _scaleX +
                                 (i * (_bounds.width() / _scaleX) / numXTicks));
            _xTicksText.push_back(addNonScalableText(QString::number(x)));
            _xTicksText.last()->setPos(x * _scaleX, _bounds.bottom() + 15);
        }
    }
    else
    {
        for (int i = 0; i <= numXTicks; ++i)
        {
            auto x =
                static_cast<int>(_bounds.left() / _scaleX +
                                 (i * (_bounds.width() / _scaleX) / numXTicks));
            QDateTime currentDate = _startDate.addSecs(x);
            _xTicksText.push_back(
                addNonScalableText(currentDate.toString("dd.MM.yyyy")));
            _xTicksText.last()->setPos(x * _scaleX, _bounds.bottom() + 15);
        }
    }

    for (int j = 0; j <= numYTicks; ++j)
    {
        float y     = _bounds.bottom() / _scaleY -
                      (j * (_bounds.height() / _scaleY) / numYTicks);
        float label = _bounds.top()   / _scaleY +
                      (j * (_bounds.height() / _scaleY) / numYTicks);
        _yTicksText.push_back(addNonScalableText(QString::number(label)));
        _yTicksText.last()->setPos(_bounds.left() - MARGIN / 2, y * _scaleY);
    }
}

/// Plots the graph.
void DiagramScene::drawGraph(DiagramList* list)
{
    QPainterPath path;

    if (list->getPath(path, _scaleX, _scaleY))
    {
        QPen pen(list->getColor(), 2, Qt::SolidLine, Qt::SquareCap, Qt::RoundJoin);
        pen.setCosmetic(true);
        _graphs.push_back(addPath(path, pen));
        addCaption(list->getName(), pen);

        int last = _graphs.size() - 1;

        /**
         * For correct display the graph needs to be flipped vertically and then
         * translated back to its original position
         */
        int verticalShift =
                static_cast<int>(2 *
                                 (list->minYValue() *
                                  _scaleY) + (_graphs[last]->boundingRect()).height());
        _graphs[last]->setTransform(QTransform(QMatrix(1,0,0,-1,0,verticalShift)));
    }
}

/// Returns the y-value at which the x-axis should cross the y-axis.
/// This value is zero if minYValue<0<maxYValue and minYValue otherwise.
int DiagramScene::getXAxisOffset()
{
    return (_bounds.top() <= 0 && _bounds.bottom() > 0) ? (int)(_bounds.bottom() +
                                                                _bounds.top()) : (int)_bounds.
           bottom();
}

/// Returns the x-value at which the y-axis should cross the x-axis.
/// This value is zero if minXValue<0<maxXValue and minXValue otherwise.
int DiagramScene::getYAxisOffset()
{
    return (_bounds.left() <= 0 && _bounds.right() > 0) ? 0 : (int)_bounds.left();
}

/// Initialises the coordinate axes, adds labels and/or units to the axes,
/// calculates axes-lengths, offsets, etc.
void DiagramScene::initialize()
{
    QPen pen(Qt::black, 2, Qt::SolidLine, Qt::SquareCap, Qt::RoundJoin);
    pen.setCosmetic(true);

    setXAxis(addArrow(_bounds.width(),  0, pen));
    setYAxis(addArrow(_bounds.height(), -90, pen));
    _xLabel = addNonScalableText(" ");
    _yLabel = addNonScalableText(" ");
    _yLabel->setRotation(-90);

    _xUnit = addNonScalableText(" ");
    _yUnit = addNonScalableText(" ");

    update();
}

/// Updates the (unscaled) boundaries of the visible coordinate system when a new
/// list is added (boundaries are rescaled in the constructGrid-method
void DiagramScene::setDiagramBoundaries(DiagramList* list)
{
    if (!_lists.isEmpty())
    {
        if (list->minXValue() < _unscaledBounds.left())
            _unscaledBounds.setLeft(list->minXValue());
        if (list->minYValue() < _unscaledBounds.top())
            _unscaledBounds.setTop(list->minYValue());
        if (list->maxXValue() > _unscaledBounds.right())
            _unscaledBounds.setRight(list->maxXValue());
        if (list->maxYValue() > _unscaledBounds.bottom())
            _unscaledBounds.setBottom(list->maxYValue());
        if (_startDate > list->getStartDate())
            _startDate = list->getStartDate();
    }
    else
    {
        _unscaledBounds.setRect(list->minXValue(), list->minYValue(),
                                list->maxXValue() - list->minXValue(),
                                list->maxYValue() - list->minYValue());
        _startDate = list->getStartDate();
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
    _xAxis->setPos(_bounds.left(),getXAxisOffset());
    _yAxis->setPos(getYAxisOffset(),_bounds.bottom());
    _xAxis->setLength(_bounds.width());
    _yAxis->setLength(_bounds.height());

    _xLabel->setPos( _bounds.left() + _bounds.width() / 2, _bounds.bottom() + 1.5 * MARGIN );
    _yLabel->setPos( _bounds.left() - 1.5 * MARGIN, _bounds.top() + _bounds.height() / 2 );

    _xUnit->setPos( _bounds.right(), _bounds.bottom() + 1.2 * MARGIN);
    _yUnit->setPos( _bounds.left(), _bounds.top() - 0.5 * MARGIN);

    /* update graphs and their captions */
    QRectF rect;
    for (int i = 0; i < _graphs.size(); i++)
    {
        rect = _graphs[i]->boundingRect();
        auto offset = static_cast<int>(fabs(rect.bottom() - _bounds.bottom()) -
                                       fabs(rect.top() - _bounds.top()));
        _graphs[i]->setPos(0, offset);

        rect = itemsBoundingRect();
        _graphCaptions[i]->setPos(_bounds.left(),rect.bottom() + 10);
    }
}
