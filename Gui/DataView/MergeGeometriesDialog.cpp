/**
 * \file
 * \author Karsten Rink
 * \date   2013-05-29
 * \brief  Implementation of the MergeGeometriesDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MergeGeometriesDialog.h"

#include "GEOObjects.h"
#include <QStringList>
#include <QStringListModel>

#include "OGSError.h"

MergeGeometriesDialog::MergeGeometriesDialog(GeoLib::GEOObjects* geoObjects, QDialog* parent)
	: QDialog(parent), _geo_objects(geoObjects), _allGeo(new QStringListModel), _selGeo(new QStringListModel)
{
	setupUi(this);

	std::vector<std::string> geoNames;
	_geo_objects->getGeometryNames(geoNames);

	// get station names
	std::vector<std::string> geo_station_names;
	_geo_objects->getStationVectorNames(geo_station_names);

	// merge method does currently not merge stations, converter function needed first
	//geoNames.reserve(geo_station_names.size());
	//std::copy(geo_station_names.begin(), geo_station_names.end(), std::back_inserter(geoNames));

	size_t nGeoObjects(geoNames.size());

	QStringList list;
	for (unsigned i = 0; i < nGeoObjects; ++i)
		list.append(QString::fromStdString(geoNames[i]));

	if (list.empty())
	{
		this->selectGeoButton->setDisabled(true);
		this->deselectGeoButton->setDisabled(true);
		list.append("(No geometry available.)");
	}
	_allGeo->setStringList(list);
	this->allGeoView->setModel(_allGeo);
	this->selectedGeoView->setModel(_selGeo);

	std::string new_geo_name("MergedGeometry");
	_geo_objects->isUniquePointVecName(new_geo_name);
	this->newGeoNameEdit->setText(QString::fromStdString(new_geo_name));
}

MergeGeometriesDialog::~MergeGeometriesDialog()
{
	delete _allGeo;
	delete _selGeo;
}

void MergeGeometriesDialog::on_selectGeoButton_pressed()
{
	QModelIndexList selected = this->allGeoView->selectionModel()->selectedIndexes();
	QStringList list = _selGeo->stringList();

	for (QModelIndexList::iterator it = selected.begin(); it != selected.end(); ++it)
	{
		list.append(it->data().toString());

		_allGeo->removeRow(it->row());
	}
	_selGeo->setStringList(list);
}

void MergeGeometriesDialog::on_deselectGeoButton_pressed()
{
	QModelIndexList selected = this->selectedGeoView->selectionModel()->selectedIndexes();
	QStringList list = _allGeo->stringList();

	for (QModelIndexList::iterator it = selected.begin(); it != selected.end(); ++it)
	{
		list.append(it->data().toString());

		_selGeo->removeRow(it->row());
	}
	_allGeo->setStringList(list);
}


void MergeGeometriesDialog::accept()
{
	std::vector<std::string> selected_geo_objects = this->getSelectedGeometries(_selGeo->stringList());
	std::string new_geometry_name = this->newGeoNameEdit->text().toStdString();
	
	if (selected_geo_objects.size()>1)
	{
		int result = _geo_objects->mergeGeometries(selected_geo_objects, new_geometry_name);
		if (result<0) 
			OGSError::box("Points are missing for\n at least one geometry.");
		else
			this->done(QDialog::Accepted);
	}
	else
		OGSError::box("At least two geometries need\n to be selected for merging.");
}

void MergeGeometriesDialog::reject()
{
	this->done(QDialog::Rejected);
}

std::vector<std::string> MergeGeometriesDialog::getSelectedGeometries(QStringList list)
{
	std::vector<std::string> indexList;
	for (QStringList::iterator it = list.begin(); it != list.end(); ++it)
		indexList.push_back(it->toStdString());
	return indexList;
}
