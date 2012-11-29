/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file ListPropertiesDialog.h
 *
 * Created on by Karsten Rink
 */

#ifndef LISTPROPERTIESDIALOG_H
#define LISTPROPERTIESDIALOG_H

#include "GEOModels.h"
#include "Station.h"
#include <QDialog>
#include <vector>

class QLabel;
class QLineEdit;
class QDialogButtonBox;

/**
 * \brief A dialog for selecting a subset of a station list based on the properties of that list.
 *
 * A dialog for selecting a subset of a station list based on the properties of that list.
 * Note: Currently, this dialog only works if the list is loaded from a database.
 */
class ListPropertiesDialog : public QDialog
{
	Q_OBJECT

public:
	ListPropertiesDialog(std::string listName, GEOModels* geoModels, QDialog* parent = 0);
	~ListPropertiesDialog();

private:
	int getPropertyBounds(const std::vector<GeoLib::Point*>* stations,
	                      const std::string &prop,
	                      double &minVal,
	                      double &maxVal);
	void setupDialog();

	QDialogButtonBox* _buttonBox;   /// The buttons used in this dialog.
	std::vector<QLabel*> _propLabel; /// The names of the properties.
	std::vector<QLineEdit*> _minValue; /// The minimum values of each property.
	std::vector<QLineEdit*> _maxValue; /// The maximum values of each property.

	std::string _listName;
	GEOModels* _geoModels;

private slots:
	void accept();
	void reject();

signals:
	void propertyBoundariesChanged(std::string name, std::vector<PropertyBounds> bounds);
};

#endif //LISTPROPERTIESDIALOG_H
