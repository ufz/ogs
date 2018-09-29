/**
 * \file
 * \author Karsten Rink
 * \date   2010-03-16
 * \brief  Implementation of the StratWindow class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Station.h"
#include "StratWindow.h"

StratWindow::StratWindow(GeoLib::StationBorehole* station,
                         std::map<std::string, DataHolderLib::Color>* stratColors,
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
