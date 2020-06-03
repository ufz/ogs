/**
 * \file
 * \author Karsten Rink
 * \date   2010-01-21
 * \brief  Implementation of the GMSHPrefsDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// Base
#include "StringTools.h"

// Qt/Base
#include "StrictDoubleValidator.h"
#include "StrictIntValidator.h"

#include "GEOObjects.h"
#include "GMSHPrefsDialog.h"
#include <QStringList>
#include <QStringListModel>

#include "OGSError.h"

GMSHPrefsDialog::GMSHPrefsDialog(GeoLib::GEOObjects const& geoObjects, QDialog* parent)
    : QDialog(parent), allGeo_(new QStringListModel), selGeo_(new QStringListModel)
{
    setupUi(this);

    // default parameters
    this->param1->setText("2");
    this->param2->setText("0.3");
    this->param3->setText("0.05");
    this->param4->setText("0");

    // object will be deleted by Qt
    auto* max_number_of_points_in_quadtree_leaf_validator(
        new StrictIntValidator(1, 1000, this->param1));
    param1->setValidator (max_number_of_points_in_quadtree_leaf_validator);
    // object will be deleted by Qt
    auto* mesh_density_scaling_pnts_validator(
        new StrictDoubleValidator(0, 1, 5, this->param2));
    param2->setValidator (mesh_density_scaling_pnts_validator);
    // object will be deleted by Qt#
    auto* mesh_density_scaling_stations_validator(
        new StrictDoubleValidator(0, 1, 5, this->param3));
    param3->setValidator (mesh_density_scaling_stations_validator);

    std::vector<std::string> geoNames;
    geoObjects.getGeometryNames(geoNames);

    // get station names
    std::vector<std::string> geo_station_names;
    geoObjects.getStationVectorNames(geo_station_names);

    std::copy(geo_station_names.begin(), geo_station_names.end(),
              geoNames.begin());

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
        list.append("[No geometry available.]");
    }
    allGeo_->setStringList(list);
    this->allGeoView->setModel(allGeo_);
    this->selectedGeoView->setModel(selGeo_);
    this->radioAdaptive->toggle(); // default is adaptive meshing
    this->on_radioAdaptive_toggled(true);
}

GMSHPrefsDialog::~GMSHPrefsDialog()
{
    delete allGeo_;
    delete selGeo_;
}

void GMSHPrefsDialog::on_selectGeoButton_pressed()
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

void GMSHPrefsDialog::on_deselectGeoButton_pressed()
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

void GMSHPrefsDialog::on_radioAdaptive_toggled(bool isTrue)
{
    if (isTrue) // meshing set to adaptive
    {
        this->param1->setEnabled(true);
        this->param2->setEnabled(true);
        this->param3->setEnabled(true);
        this->param4->setEnabled(false);
    }
    else // meshing set to homogeneous
    {
        this->param1->setEnabled(false);
        this->param2->setEnabled(false);
        this->param3->setEnabled(false);
        this->param4->setEnabled(true);
    }
}

void GMSHPrefsDialog::accept()
{
    if (this->selGeo_->stringList().empty())
    {
        OGSError::box("No geometry selected. Geometric data\n is necessary for mesh generation.");
        return;
    }

    std::vector<std::string> selectedObjects = this->getSelectedObjects(selGeo_->stringList());
    unsigned max_number_of_points_in_quadtree_leaf (10);
    double mesh_density_scaling_pnts(0.5);
    double mesh_density_scaling_stations (0.05);
    double val4(-1);

    if (this->radioAdaptive->isChecked())
    {
        double const min_scaling_factor (1e-10);
        max_number_of_points_in_quadtree_leaf =
            BaseLib::str2number<unsigned>(param1->text().toStdString());
        if (max_number_of_points_in_quadtree_leaf == 0)
        {
            max_number_of_points_in_quadtree_leaf = 10;
        }
        mesh_density_scaling_pnts = fabs (param2->text().toDouble());
        if (mesh_density_scaling_pnts < min_scaling_factor)
        {
            mesh_density_scaling_pnts = min_scaling_factor;
        }
        mesh_density_scaling_stations = param3->text().toDouble();
        if (mesh_density_scaling_stations < min_scaling_factor)
        {
            mesh_density_scaling_stations = min_scaling_factor;
        }
    }
    else
    {
        val4 = param4->text().toDouble();
    }

    bool delete_geo_file = this->geoFileDelete->isChecked();
    emit requestMeshing(selectedObjects,
                        max_number_of_points_in_quadtree_leaf,
                        mesh_density_scaling_pnts,
                        mesh_density_scaling_stations,
                        val4,
                        delete_geo_file);
    this->done(QDialog::Accepted);
}

void GMSHPrefsDialog::reject()
{
    this->done(QDialog::Rejected);
}

std::vector<std::string> GMSHPrefsDialog::getSelectedObjects(QStringList list)
{
    std::vector<std::string> indexList;
    std::transform(list.begin(), list.end(), std::back_inserter(indexList),
                   [](auto const& index) { return index.toStdString(); });
    return indexList;
}
