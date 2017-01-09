/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Implementation of the DiagramView class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DiagramView.h"
#include <QGraphicsTextItem>
#include <math.h>

DiagramView::DiagramView(QWidget* parent) : QGraphicsView(parent)
{
    _scene = new DiagramScene();
    setScene(_scene);
    initialize();
}

DiagramView::DiagramView(DiagramList* list, QWidget* parent) : QGraphicsView(parent)
{
    _scene = new DiagramScene(list);
    setScene(_scene);
    initialize();
}

DiagramView::~DiagramView()
{
    delete _scene;
}

void DiagramView::addGraph(DiagramList* list)
{
    _scene->addGraph(list);
    update();
}

int DiagramView::getHeight()
{
    return static_cast<int>((_scene->itemsBoundingRect()).height());
}

int DiagramView::getWidth()
{
    return static_cast<int>((_scene->itemsBoundingRect()).width());
}

/**
 * Initialises the view.
 */
void DiagramView::initialize()
{
    //QMatrix currentMatrix = matrix();
    //setMatrix(currentMatrix * scene->getTransformationMatrix());

    setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

    update();
}

/*
 * Keeps the aspect ration of the labels when the view is resized.
 * It is only necessary to call this if
 *        Qt::AspectRatioMode == Qt::IgnoreAspectRatio.
 * Also, this method is kind of annoying because you have to set the
 * appropriate transform for every single QGraphicsTextItem separately.
 */
/*
   void DiagramView::keepItemAspectRatio()
   {
    double xFactor = transform().mapRect(QRectF(0, 0, 1, 1)).width();
    double yFactor = transform().mapRect(QRectF(0, 0, 1, 1)).height();
    QMatrix invertedScaling;
    invertedScaling.scale(1.0 , xFactor / yFactor);

    scene->xLabel->setTransform(QTransform(invertedScaling));
    scene->yLabel->setTransform(QTransform(invertedScaling));
    scene->yLabel->rotate(-90);
   }
 */

QSize DiagramView::minimumSizeHint() const
{
    return QSize(3 * _scene->MARGIN,2 * _scene->MARGIN);
}

QSize DiagramView::sizeHint() const
{
    return QSize(6 * _scene->MARGIN, 4 * _scene->MARGIN);
}

void DiagramView::resizeEvent(QResizeEvent* event)
{
    Q_UNUSED (event)
    update();
    //keepItemAspectRatio();
}

/**
 * Updates the view automatically when a new list is added or when
 * the window containing the view is resized or changes its state.
 * Basically, the methods makes sure that everything keeps looking
 * as it is supposed to.
 */
void DiagramView::update()
{
    //setResizeAnchor(QGraphicsView::AnchorViewCenter);
    QRectF viewRect = _scene->itemsBoundingRect();
    _scene->setSceneRect(viewRect);
    QRectF sceneInView(0 /*_scene->MARGIN*/,_scene->MARGIN / 2,
                       viewRect.width() /*+_scene->MARGIN*/,viewRect.height() + _scene->MARGIN);
    fitInView(sceneInView, Qt::IgnoreAspectRatio);
}
