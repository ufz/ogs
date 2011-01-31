/**
 * \file StratWindow.cpp
 * 2010/03/16 - KR Initial implementation
 */

#include "StratWindow.h"
#include "Station.h"


StratWindow::StratWindow(GEOLIB::StationBorehole* station, std::map<std::string, GEOLIB::Color*> *stratColors, QWidget* parent) : QWidget(parent)
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
	int width = (stationView->getWidth()>800) ? 800 : stationView->getWidth();
	int height = (stationView->getHeight()>600) ? 600 : stationView->getHeight();
	resize(width, height);
}
