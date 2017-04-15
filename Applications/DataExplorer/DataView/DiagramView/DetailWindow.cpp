/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Implementation of the DetailWindow class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DetailWindow.h"

#include <QFileDialog>
#include <QSettings>

#include "Applications/DataHolderLib/Color.h"
#include "DiagramPrefsDialog.h"

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

    for (auto& list : lists)
        stationView->addGraph(list);

    resizeWindow();
}

DetailWindow::DetailWindow(DiagramList* list, QWidget* parent) : QWidget(parent)
{
    setupUi(this);
    stationView->setRenderHints( QPainter::Antialiasing );
    stationView->addGraph(list);
    resizeWindow();
}

DetailWindow::DetailWindow(std::vector<std::size_t> data, QWidget* parent) : QWidget(parent)
{
    setupUi(this);
    std::size_t nEntries = data.size();
    std::vector< std::pair<float, float> > list_data(nEntries);

    for (std::size_t i=0; i<nEntries; i++)
        list_data.emplace_back(static_cast<float>(i),
                               static_cast<float>(data[i]));

    auto* list = new DiagramList();
    list->setList(list_data);
    list->setXUnit("Value");
    list->setYUnit("Amount");
    list->setName("Histogram");
    stationView->setRenderHints( QPainter::Antialiasing );
    stationView->addGraph(list);
    resizeWindow();
}

DetailWindow::~DetailWindow() = default;

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
    DataHolderLib::Color const c = DataHolderLib::getRandomColor();
    QColor colour(c[0], c[1], c[2]);
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
    QSettings settings;
    QString fileName = QFileDialog::getOpenFileName( this, "Select data file to open",
                                                     settings.value("lastOpenedFileDirectory").toString(),
                                                     "Text files (*.txt);;All files (* *.*)");
    if (!fileName.isEmpty())
    {
        QDir dir = QDir(fileName);
        settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
        auto* prefs = new DiagramPrefsDialog(fileName, this);
        prefs->setAttribute(Qt::WA_DeleteOnClose);
        prefs->show();
    }
}

