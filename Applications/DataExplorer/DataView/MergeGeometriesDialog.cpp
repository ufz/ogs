/**
 * \file
 * \author Karsten Rink
 * \date   2013-05-29
 * \brief  Implementation of the MergeGeometriesDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
    : QDialog(parent), geo_objects_(geoObjects), allGeo_(new QStringListModel), selGeo_(new QStringListModel)
{
    setupUi(this);

    std::vector<std::string> geoNames;
    geo_objects_.getGeometryNames(geoNames);

    // get station names
    std::vector<std::string> geo_station_names;
    geo_objects_.getStationVectorNames(geo_station_names);

    // merge method does currently not merge stations, converter function needed first
    //geoNames.reserve(geo_station_names.size());
    //std::copy(geo_station_names.begin(), geo_station_names.end(), std::back_inserter(geoNames));

    std::size_t nGeoObjects(geoNames.size());

    QStringList list;
    for (unsigned i = 0; i < nGeoObjects; ++i)
    {
        list.append(QString::fromStdString(geoNames[i]));
    }

    if (list.empty())
    {
        this->selectGeoButton->setDisabled(true);
        this->deselectGeoButton->setDisabled(true);
        list.append("(No geometry available.)");
    }
    allGeo_->setStringList(list);
    this->allGeoView->setModel(allGeo_);
    this->selectedGeoView->setModel(selGeo_);

    std::string new_geo_name("MergedGeometry");
    geo_objects_.isUniquePointVecName(new_geo_name);
    this->newGeoNameEdit->setText(QString::fromStdString(new_geo_name));
}

MergeGeometriesDialog::~MergeGeometriesDialog()
{
    delete allGeo_;
    delete selGeo_;
}

void MergeGeometriesDialog::on_selectGeoButton_pressed()
{
    QModelIndexList selected = this->allGeoView->selectionModel()->selectedIndexes();
    QStringList list = selGeo_->stringList();

    for (auto& index : selected)
    {
        list.append(index.data().toString());

        allGeo_->removeRow(index.row());
    }
    selGeo_->setStringList(list);
}

void MergeGeometriesDialog::on_deselectGeoButton_pressed()
{
    QModelIndexList selected = this->selectedGeoView->selectionModel()->selectedIndexes();
    QStringList list = allGeo_->stringList();

    for (auto& index : selected)
    {
        list.append(index.data().toString());

        selGeo_->removeRow(index.row());
    }
    allGeo_->setStringList(list);
}

void MergeGeometriesDialog::accept()
{
    if (selGeo_->stringList().size() > 1)
    {
        this->done(QDialog::Accepted);
    }
    else
    {
        OGSError::box(
            "At least two geometries need\n to be selected for merging.");
    }
}

void MergeGeometriesDialog::reject()
{
    this->done(QDialog::Rejected);
}

std::vector<std::string> MergeGeometriesDialog::getSelectedGeometries() const
{
    std::vector<std::string> indexList;
    QStringList const& list(selGeo_->stringList());
    std::transform(list.begin(), list.end(), std::back_inserter(indexList),
                   [](auto const& index) { return index.toStdString(); });
    return indexList;
}

std::string MergeGeometriesDialog::getGeometryName() const
{
    return this->newGeoNameEdit->text().toStdString();
}
