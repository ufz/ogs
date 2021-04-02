/**
 * \file
 * \author Karsten Rink
 * \date   2010-03-16
 * \brief  Implementation of the StratView class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "StratView.h"

#include <math.h>

#include "Station.h"

StratView::~StratView()
{
    delete _scene;
}

void StratView::setStation(
    GeoLib::StationBorehole* station,
    std::map<std::string, DataHolderLib::Color>* stratColors)
{
    _scene = new StratScene(station, stratColors);
    setScene(_scene);
    initialize();
}

void StratView::initialize()
{
    setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

    update();
}

void StratView::resizeEvent(QResizeEvent* event)
{
    Q_UNUSED(event)
    update();
}

void StratView::update()
{
    QRectF viewRect = _scene->itemsBoundingRect();
    _scene->setSceneRect(viewRect);
    QRectF sceneInView(StratScene::MARGIN, StratScene::MARGIN,
                       viewRect.width() + 2 * StratScene::MARGIN,
                       viewRect.height() + 2 * StratScene::MARGIN);
    fitInView(sceneInView, Qt::IgnoreAspectRatio);
}

void StratView::saveAsImage(QString fileName)
{
    this->update();

    QRectF viewRect = _scene->itemsBoundingRect();
    QImage img(static_cast<int>(viewRect.width()) + 2 * StratScene::MARGIN, 600,
               QImage::Format_ARGB32);
    QPainter painter(&img);

    this->render(&painter);
    img.save(fileName);
}
