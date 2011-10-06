/**
 * \file DetailWindow.cpp
 * KR Initial implementation
 */

#include "Color.h"
#include "DetailWindow.h"
#include "DiagramPrefsDialog.h"

#include <QFileDialog>
#include <QSettings>

DetailWindow::DetailWindow(QWidget* parent) : QWidget(parent)
{
	setupUi(this);
	stationView->setRenderHints( QPainter::Antialiasing );

/*
    DiagramList* list  = new DiagramList();
    DiagramList* list2 = new DiagramList();


    // ==================================================
    // input files should be defined in WaterML
    // inserting the details below into the list-objects
    // kind of simulates the information that would be
    // included in a WaterML-file and is needed for
    // display
    // ==================================================

    // make up list-object for the first test station
    list->setName("Water Level Observation Station: Halberstadt 2002");
    list->setXLabel("Time");
    list->setYLabel("Water Level");
    list->setXUnit("day");
    list->setYUnit("metres");
    list->setColor(QColor(Qt::red));
    list->readList("c:\\project\\timeseries-a.stn");

    // make up list-object for the second test station
    list2->setName("Water Level Observation Station: Oschersleben 2002");
    list2->setXLabel("Time");
    list2->setYLabel("Water Level");
    list2->setXUnit("day");
    list2->setYUnit("metres");
    list2->setColor(QColor(Qt::green));
    list2->readList("c:\\project\\timeseries-b.stn");

    // ==================================================


    stationView->addGraph(list);
    stationView->addGraph(list2);

    resizeWindow();
 */
}

DetailWindow::DetailWindow(QString filename, QWidget* parent) : QWidget(parent)
{
	setupUi(this);
	stationView->setRenderHints( QPainter::Antialiasing );

	std::vector<DiagramList*> lists;
	DiagramList::readList(filename, lists);

	for (size_t i = 0; i < lists.size(); i++)
		stationView->addGraph(lists[i]);

	resizeWindow();
}

DetailWindow::DetailWindow(DiagramList* list, QWidget* parent) : QWidget(parent)
{
	setupUi(this);
	stationView->setRenderHints( QPainter::Antialiasing );
	stationView->addGraph(list);
	resizeWindow();
}

DetailWindow::~DetailWindow()
{
}

void DetailWindow::on_closeButton_clicked()
{
	this->close();
}

void DetailWindow::resizeWindow()
{
	int width = (stationView->getWidth() > 800) ? 800 : stationView->getWidth();
	int height = (stationView->getHeight() > 600) ? 600 : stationView->getHeight();
	resize(width, height);
}

void DetailWindow::addList(DiagramList* list)
{
	GEOLIB::Color* c = GEOLIB::getRandomColor();
	QColor colour((*c)[0], (*c)[1], (*c)[2]);
	delete c;
	this->addList(list, colour);
	resizeWindow();
}

void DetailWindow::addList(DiagramList* list, QColor c)
{
	list->setColor(c);
	this->stationView->addGraph(list);
}

void DetailWindow::on_addDataButton_clicked()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName = QFileDialog::getOpenFileName( this, "Select data file to open",
	                                                 settings.value(
	                                                         "lastOpenedFileDirectory").
	                                                 toString(),
	                                                 "Text files (*.txt);;All files (* *.*)");
	if (!fileName.isEmpty())
	{
		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
		DiagramPrefsDialog* prefs = new DiagramPrefsDialog(fileName, this);
		prefs->setAttribute(Qt::WA_DeleteOnClose);
		prefs->show();
	}
}

