/**
 * \file
 * \author Karsten Rink
 * \date   2012-04-17
 * \brief  Implementation of the LinearEditDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LinearEditDialog.h"

LinearEditDialog::LinearEditDialog(const GeoLib::Polyline &line,
                                   const std::vector<std::size_t> &dis_nodes,
                                   const std::vector<double> &dis_values,
                                   QDialog* parent)
    : QDialog(parent), _line(line)
{
    setupUi(this);
    setupDialog(dis_nodes, dis_values);
}

void LinearEditDialog::setupDialog(const std::vector<std::size_t> &dis_nodes,
                                   const std::vector<double> &dis_values)
{
    std::size_t nPoints(_line.getNumberOfPoints());
    this->tableWidget->setRowCount(nPoints);
    QList<QString> indexlist;

    for (std::size_t i = 0; i < nPoints; i++)
    {
        indexlist.push_back(QString::number(i));
        QTableWidgetItem *newItem = new QTableWidgetItem("");
        tableWidget->setItem(i, 0, newItem);
    }
    QStringList vHeaders(indexlist);
    tableWidget->setVerticalHeaderLabels(vHeaders);

    std::size_t nValues (dis_values.size());
    for (std::size_t i = 0; i < nValues; i++)
        tableWidget->item(dis_nodes[i],0)->setText(QString::number(dis_values[i]));
}

LinearEditDialog::~LinearEditDialog() = default;

void LinearEditDialog::on_comboBox_currentIndexChanged(int index)
{
    if (index > 0) //elevation
    {
        std::size_t nRows = tableWidget->rowCount();
        for (std::size_t i = 0; i < nRows; i++)
            tableWidget->item(i,0)->setText(QString::number(_line.getPoint(i)->getCoords()[2]));
    }
}

void LinearEditDialog::accept()
{
    std::vector< std::pair<std::size_t,double> > linear_values;

    std::size_t nRows = tableWidget->rowCount();
    for (std::size_t i = 0; i < nRows; i++)
    {
        QString row_text (tableWidget->item(i,0)->text());
        if (row_text.length() > 0)
            linear_values.push_back( std::pair<std::size_t, double>(i, row_text.toDouble()) );
    }

    emit transmitDisValues(linear_values);
    this->done(QDialog::Accepted);
}

void LinearEditDialog::reject()
{
    this->done(QDialog::Rejected);
}
