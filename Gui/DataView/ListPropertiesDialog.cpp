/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Implementation of the ListPropertiesDialog class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DateTools.h"
#include "ListPropertiesDialog.h"
#include "PropertyBounds.h"
#include "StringTools.h"
#include <QDialogButtonBox>
#include <QGridLayout>
#include <QLabel>
#include <QLineEdit>
#include <QLocale>

/**
 * Creates a new dialog.
 */
ListPropertiesDialog::ListPropertiesDialog(std::string listName,
                                           GEOModels* geoModels,
                                           QDialog* parent) :
	QDialog(parent), _listName(listName), _geoModels(geoModels)
{
	setupDialog();
	show();
}

ListPropertiesDialog::~ListPropertiesDialog()
{
	for (size_t i = 0; i < _propLabel.size(); i++)
	{
		delete _propLabel[i];
		delete _minValue[i];
		delete _maxValue[i];
	}

	delete _buttonBox;
}

/// Constructs a dialog window based on the properties retrieved from the station objects
void ListPropertiesDialog::setupDialog()
{
	int i = 0;
	double minVal = 0, maxVal = 0;

	const std::vector<GeoLib::Point*>* stations ( _geoModels->getStationVec(_listName));

	std::map<std::string,
	         double> properties = static_cast<GeoLib::Station*>((*stations)[0])->getProperties();
	QGridLayout* layout = new QGridLayout;

	setWindowTitle("List Properties");

	for(std::map<std::string, double>::const_iterator it = properties.begin();
	    it != properties.end(); ++it)
	{
		QLabel* _prop = new QLabel(this);
		QLineEdit* _min = new QLineEdit(this);
		QLineEdit* _max = new QLineEdit(this);
		_prop->setText(QString::fromStdString(it->first));

		if (getPropertyBounds(stations, it->first, minVal, maxVal))
		{
			_min->setText(QString::number(minVal, 'f'));
			if (_prop->text().compare("date") == 0)
				_min->setText(QString::fromStdString(BaseLib::date2string(minVal)));

			_max->setText(QString::number(maxVal, 'f'));
			if (_prop->text().compare("date") == 0)
				_max->setText(QString::fromStdString(BaseLib::date2string(maxVal)));
		}

		_propLabel.push_back(_prop);
		_minValue.push_back(_min);
		_maxValue.push_back(_max);

		layout->addWidget( _propLabel[i], i, 0 );
		layout->addWidget( _minValue[i], i, 1 );
		layout->addWidget( _maxValue[i], i, 2 );
		i++;
	}

	_buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
	connect(_buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
	connect(_buttonBox, SIGNAL(rejected()), this, SLOT(reject()));

	layout->addWidget(_buttonBox, i + 1, 1, 1, 2 );

	setLayout(layout);
}

int ListPropertiesDialog::getPropertyBounds(const std::vector<GeoLib::Point*>* stations,
                                            const std::string &prop,
                                            double &minVal,
                                            double &maxVal)
{
	if (!stations->empty())
	{
		std::map<std::string, double> properties (
		        static_cast<GeoLib::Station*>((*stations)[0])->getProperties());
		minVal = properties[prop];
		maxVal = properties[prop];

		size_t size = stations->size();
		for (size_t i = 1; i < size; i++)
		{
			properties = static_cast<GeoLib::Station*>((*stations)[i])->getProperties();
			if (minVal > properties[prop])
				minVal = properties[prop];
			if (maxVal < properties[prop])
				maxVal = properties[prop];
		}
		return 1;
	}
	return 0;
}

/// Instructions if the OK-Button has been pressed.
void ListPropertiesDialog::accept()
{
	std::vector<PropertyBounds> bounds;
	int noProp = _propLabel.size();
	double minVal, maxVal;

	for (int i = 0; i < noProp; i++)
	{
		if (_propLabel[i]->text().compare("date") == 0)
		{
			minVal = BaseLib::xmlDate2int(_minValue[i]->text().toStdString());
			maxVal = BaseLib::xmlDate2int(_maxValue[i]->text().toStdString());
		}
		else
		{
			minVal = _minValue[i]->text().toDouble();
			maxVal = _maxValue[i]->text().toDouble();
		}
		PropertyBounds b(_propLabel[i]->text().toStdString(), minVal, maxVal);
		bounds.push_back(b);
	}

	emit propertyBoundariesChanged(_listName, bounds);

	this->done(QDialog::Accepted);
}

/// Instructions if the Cancel-Button has been pressed.
void ListPropertiesDialog::reject()
{
	this->done(QDialog::Rejected);
}
