/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file LinearEditDialog.cpp
 *
 * Created on 2012-04-17 by Karsten Rink
 */

#include "LinearEditDialog.h"

LinearEditDialog::LinearEditDialog(const GeoLib::Polyline &line, const std::vector<size_t> &dis_nodes, const std::vector<double> &dis_values, QDialog* parent)
	: QDialog(parent), _line(line)
{
	setupUi(this);
	setupDialog(dis_nodes, dis_values);
}

void LinearEditDialog::setupDialog(const std::vector<size_t> &dis_nodes, const std::vector<double> &dis_values)
{
	size_t nPoints(_line.getNumberOfPoints());
	this->tableWidget->setRowCount(nPoints);
	QList<QString> indexlist;

	for (size_t i=0; i<nPoints; i++)
	{
		indexlist.push_back(QString::number(i));
		QTableWidgetItem *newItem = new QTableWidgetItem("");
		tableWidget->setItem(i, 0, newItem);
	}
	QStringList vHeaders(indexlist);
	tableWidget->setVerticalHeaderLabels(vHeaders);

	size_t nValues (dis_values.size());
	for (size_t i=0; i<nValues; i++)
		tableWidget->item(dis_nodes[i],0)->setText(QString::number(dis_values[i]));
}

LinearEditDialog::~LinearEditDialog()
{
}

void LinearEditDialog::on_comboBox_currentIndexChanged(int index)
{
	if (index>0) //elevation
	{
		size_t nRows = tableWidget->rowCount();
		for (size_t i=0; i<nRows; i++)
			tableWidget->item(i,0)->setText(QString::number(_line[i]->getCoords()[2]));
	}
}

void LinearEditDialog::accept()
{
	std::vector< std::pair<size_t,double> > linear_values;

	size_t nRows = tableWidget->rowCount();
	for (size_t i=0; i<nRows; i++)
	{
		QString row_text (tableWidget->item(i,0)->text());
		if (row_text.length()>0)
			linear_values.push_back( std::pair<size_t, double>(i, row_text.toDouble()) );
	}

	emit transmitDisValues(linear_values);
	this->done(QDialog::Accepted);
}

void LinearEditDialog::reject()
{
	this->done(QDialog::Rejected);
}
