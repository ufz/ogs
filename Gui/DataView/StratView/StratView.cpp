/**
 * \file StratView.cpp
 * 2010/03/16 - KR Initial implementation
 */

#include "Station.h"
#include "StratView.h"
#include <math.h>

StratView::~StratView()
{
	delete _scene;
}

void StratView::setStation(GEOLIB::StationBorehole* station,
                           std::map<std::string, GEOLIB::Color*>* stratColors)
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
	Q_UNUSED (event)
	update();
}

void StratView::update()
{
	QRectF viewRect = _scene->itemsBoundingRect();
	_scene->setSceneRect(viewRect);
	QRectF sceneInView(_scene->MARGIN,_scene->MARGIN,
	                   viewRect.width() + 2 * _scene->MARGIN,
	                   viewRect.height() + 2 * _scene->MARGIN);
	fitInView(sceneInView, Qt::IgnoreAspectRatio);
}

void StratView::saveAsImage(QString fileName)
{
	this->update();

	QRectF viewRect = _scene->itemsBoundingRect();
	QImage img(
	        static_cast<int>(viewRect.width()) + 2 * _scene->MARGIN, 600, QImage::Format_ARGB32);
	QPainter painter(&img);

	this->render(&painter);
	img.save(fileName);
}
