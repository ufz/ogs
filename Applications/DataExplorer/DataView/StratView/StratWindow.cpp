// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "DataView/StratView/StratWindow.h"

#include "GeoLib/Station.h"

StratWindow::StratWindow(
    GeoLib::StationBorehole* station,
    std::map<std::string, DataHolderLib::Color>* stratColors,
    QWidget* parent)
    : QWidget(parent)
{
    setupUi(this);
    stationView->setRenderHints(QPainter::Antialiasing);
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
    int height =
        (stationView->getHeight() > 600) ? 600 : stationView->getHeight();
    resize(width, height);
}
