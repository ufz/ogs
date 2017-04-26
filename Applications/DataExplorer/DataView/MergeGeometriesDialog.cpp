/**
 * \file
 * \author Karsten Rink
 * \date   2013-05-29
 * \brief  Implementation of the MergeGeometriesDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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

MergeGeometriesDialog::MergeGeometriesDialog(GeoLib::GEOObjects& geoObjects, QDialog* parent)
    : QDialog(parent), _geo_objects(geoObjects), _allGeo(new QStringListModel), _selGeo(new QStringListModel)
{
    setupUi(this);

    std::vector<std::string> geoNames;
    _geo_objects.getGeometryNames(geoNames);

    // get station names
    std::vector<std::string> geo_station_names;
    _geo_objects.getStationVectorNames(geo_station_names);

    // merge method does currently not merge stations, converter function needed first
    //geoNames.reserve(geo_station_names.size());
    //std::copy(geo_station_names.begin(), geo_station_names.end(), std::back_inserter(geoNames));

    std::size_t nGeoObjects(geoNames.size());

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
    _geo_objects.isUniquePointVecName(new_geo_name);
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

    for (auto& index : selected)
    {
        list.append(index.data().toString());

        _allGeo->removeRow(index.row());
    }
    _selGeo->setStringList(list);
}

void MergeGeometriesDialog::on_deselectGeoButton_pressed()
{
    QModelIndexList selected = this->selectedGeoView->selectionModel()->selectedIndexes();
    QStringList list = _allGeo->stringList();

    for (auto& index : selected)
    {
        list.append(index.data().toString());

        _selGeo->removeRow(index.row());
    }
    _allGeo->setStringList(list);
}

void MergeGeometriesDialog::accept()
{
    if (_selGeo->stringList().size() > 1)
        this->done(QDialog::Accepted);
    else
        OGSError::box("At least two geometries need\n to be selected for merging.");
}

void MergeGeometriesDialog::reject()
{
    this->done(QDialog::Rejected);
}

std::vector<std::string> const MergeGeometriesDialog::getSelectedGeometries() const
{
    std::vector<std::string> indexList;
    QStringList const& list (_selGeo->stringList());
    for (const auto& index : list)
        indexList.push_back(index.toStdString());
    return indexList;
}

std::string MergeGeometriesDialog::getGeometryName() const
{
    return this->newGeoNameEdit->text().toStdString();
}
