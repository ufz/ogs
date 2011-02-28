/**
 * \file DetailWindow.cpp
 * KR Initial implementation
 */

#include "DetailWindow.h"

/**
 * Creates an empty window.
 */
DetailWindow::DetailWindow(QWidget* parent) : QWidget(parent)
{
	setupUi(this);
	stationView->setRenderHints( QPainter::Antialiasing );

	DiagramList* list  = new DiagramList();
	DiagramList* list2 = new DiagramList();


	/* ================================================== */
    /* input files should be defined in WaterML */
	/* inserting the details below into the list-objects */
	/* kind of simulates the information that would be */
	/* included in a WaterML-file and is needed for */
	/* display */
	/* ================================================== */

	/* make up list-object for the first test station *
	list->setName("Water Level Observation Station: Halberstadt 2002" /*Oschersleben 2003*);
	list->setXLabel("Time");
	list->setYLabel("Water Level");
	list->setXUnit("day");
	list->setYUnit("metres");
	list->setColor(QColor(Qt::red));
	list->readList("c:\\project\\timeseries-a.stn");

	/* make up list-object for the second test station *
	list2->setName("Water Level Observation Station: Oschersleben 2002");
	list2->setXLabel("Time");
	list2->setYLabel("Water Level");
	list2->setXUnit("day");
	list2->setYUnit("metres");
	list2->setColor(QColor(Qt::green));
	list2->readList("c:\\project\\timeseries-b.stn");

	/* ================================================== */


	stationView->addGraph(list);
	stationView->addGraph(list2);

	resizeWindow();
}

DetailWindow::DetailWindow(QString filename, QWidget* parent) : QWidget(parent)
{
	setupUi(this);
	stationView->setRenderHints( QPainter::Antialiasing );

	DiagramList::readList(filename, _list);

	for (size_t i=0; i<_list.size(); i++)
		stationView->addGraph(_list[i]);

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
	for (size_t i=0; i<_list.size(); i++)
		delete _list[i];
}

void DetailWindow::on_closeButton_clicked()
{
	this->close();
}


void DetailWindow::resizeWindow()
{
	int width = (stationView->getWidth()>800) ? 800 : stationView->getWidth();
	int height = (stationView->getHeight()>600) ? 600 : stationView->getHeight();
	resize(width, height);
}
