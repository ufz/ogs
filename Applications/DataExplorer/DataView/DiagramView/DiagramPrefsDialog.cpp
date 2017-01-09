/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Implementation of the DiagramPrefsDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DetailWindow.h"
#include "DiagramList.h"
#include "DiagramPrefsDialog.h"
#include "OGSError.h"
#include "Station.h"

#include <QCheckBox>
#include <QFileDialog>
#include <QMessageBox>

DiagramPrefsDialog::DiagramPrefsDialog(const GeoLib::Station* stn,
                                       const QString &listName,
                                       //DatabaseConnection* db,
                                       QDialog* parent)
    : QDialog(parent), _listID(-1), _stationID(-1), _window(nullptr)
{
    setAttribute(Qt::WA_DeleteOnClose);

    setupUi(this);
    stationNameLabel->setText(QString::fromStdString(stn->getName()));
    stationTypeLabel->setText(listName);
}

DiagramPrefsDialog::DiagramPrefsDialog(GeoLib::Station* stn, QDialog* parent)
    : QDialog(parent), _listID(-1), _stationID(-1), _window(nullptr)
{
    setupUi(this);
    stationNameLabel->setText(QString::fromStdString(stn->getName()));
    stationTypeLabel->setText("");
    DiagramList::readList(stn->getSensorData(), _list);

    fromDateLine->setText(QString::number(stn->getSensorData()->getStartTime()));
    toDateLine->setText(QString::number(stn->getSensorData()->getEndTime()));
    this->createVisibilityCheckboxes();
}

DiagramPrefsDialog::DiagramPrefsDialog(const QString &filename,
                                       DetailWindow* window,
                                       QDialog* parent)
    : QDialog(parent), _listID(-1), _stationID(-1), _window(window)
{
    QFileInfo fi(filename);
    setupUi(this);
    stationNameLabel->setText(fi.baseName());
    stationTypeLabel->setText("");
    this->loadFile(filename);
}

DiagramPrefsDialog::~DiagramPrefsDialog()
{
    this->destroy();
}

void DiagramPrefsDialog::accept()
{
    if ((fromDateLine->text().length() > 0) && (toDateLine->text().length() > 0) &&
        (!_list.empty()))
    {
        // data has been loaded
        if (_list[0]->size() > 0)
        {
            bool window_is_empty(false);
            if (_window == NULL)
            {
                _window = new DetailWindow();
                _window->setAttribute(Qt::WA_DeleteOnClose);
                window_is_empty = true;
            }

            for (std::size_t i = 0; i < _list.size(); i++)
                if (this->_visability[i]->isChecked())
                {
                    _window->addList(_list[i]);
                    window_is_empty = false;
                }

            if (!window_is_empty)
            {
                _window->show();
                this->done(QDialog::Accepted);
            }
            else
            {
                delete _window;
                _window = NULL;
                OGSError::box("No dataset selected.");
            }
        }
        else
        {
            OGSError::box("Invalid station data.");
            this->done(QDialog::Rejected);
        }
    }
    else
        OGSError::box("No data found...");
}

void DiagramPrefsDialog::reject()
{
    this->done(QDialog::Rejected);
}

void DiagramPrefsDialog::on_loadFileButton_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this,
                                                    "Select time series file to open",
                                                    "",
                                                    "Time series files (*.stn *.txt)");
    if (!fileName.isEmpty())
        loadFile(fileName);
}

int DiagramPrefsDialog::loadFile(const QString &filename)
{
    if (DiagramList::readList(filename, _list))
    {
        for (std::size_t i = 0; i < _list.size(); i++)
        {
            //_list[i]->setName(stationTypeLabel->text() + ": " + stationNameLabel->text());
            _list[i]->setXLabel("Time");
            //_list[i]->setYLabel("Water Level");
            _list[i]->setXUnit("day");
            //_list[i]->setYUnit("metres");
            _list[i]->setColor(QColor(Qt::red));
        }
        fromDateLine->setText(_list[0]->getStartDate().toString("dd.MM.yyyy")); //QString::number(_list[0]->minXValue()));
        QDateTime endDate =
                _list[0]->getStartDate().addSecs(static_cast<int>(_list[0]->maxXValue()));
        toDateLine->setText(endDate.toString("dd.MM.yyyy")); //QString::number(_list[0]->maxXValue()));
        this->createVisibilityCheckboxes();
        return 1;
    }

    OGSError::box("Error reading file.");
    return 0;
}

int DiagramPrefsDialog::loadList(const std::vector< std::pair<QDateTime, float> > &coords)
{
    if (!coords.empty())
    {
        DiagramList* l = new DiagramList;
        l->setName(stationTypeLabel->text() + ": " + stationNameLabel->text());
        l->setXLabel("Time");
        //l->setYLabel("Water Level");
        l->setXUnit("day");
        //l->setYUnit("metres");
        l->setColor(QColor(Qt::red));
        l->setList(coords);
        _list.push_back(l);
        return 1;
    }
    return 0;
}

void DiagramPrefsDialog::createVisibilityCheckboxes()
{
    for (std::size_t i = 0; i < _list.size(); i++)
    {
        QCheckBox* box = new QCheckBox(_list[i]->getName());
        box->setChecked(true);
        this->CheckBoxLayout->addWidget(box);
        _visability.push_back(box);
    }
}

