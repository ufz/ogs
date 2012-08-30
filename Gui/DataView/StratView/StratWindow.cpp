/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file StratWindow.cpp
 *
 * Created on 2010-03-16 by Karsten Rink
 */

#include "Station.h"
#include "StratWindow.h"

StratWindow::StratWindow(GeoLib::StationBorehole* station,
                         std::map<std::string, GeoLib::Color*>* stratColors,
                         QWidget* parent) : QWidget(parent)
{
	setupUi(this);
	stationView->setRenderHints( QPainter::Antialiasing );
	stationView->setStation(station, stratColors);
	resizeWindow();
}

void StratWindow::on_closeButton_clicked()
{
	this->close();
}

void StratWindow::resizeWindow()
{
	int width = (stationView->getWidth() > 800) ? 800 : stationView->getWidth();
	int height = (stationView->getHeight() > 600) ? 600 : stationView->getHeight();
	resize(width, height);
}
