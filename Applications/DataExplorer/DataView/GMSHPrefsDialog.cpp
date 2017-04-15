/**
 * \file
 * \author Karsten Rink
 * \date   2010-01-21
 * \brief  Implementation of the GMSHPrefsDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
    : QDialog(parent), _allGeo(new QStringListModel), _selGeo(new QStringListModel)
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
        new StrictDoubleValidator(1e-10, 1.0, 5, this->param2));
    param2->setValidator (mesh_density_scaling_pnts_validator);
    // object will be deleted by Qt#
    auto* mesh_density_scaling_stations_validator(
        new StrictDoubleValidator(1e-10, 1.0, 5, this->param3));
    param3->setValidator (mesh_density_scaling_stations_validator);

    std::vector<std::string> geoNames;
    geoObjects.getGeometryNames(geoNames);

    // get station names
    std::vector<std::string> geo_station_names;
    geoObjects.getStationVectorNames(geo_station_names);

    for (unsigned k(0); k < geo_station_names.size(); ++k)
        geoNames.push_back (geo_station_names[k]);

    std::size_t nGeoObjects(geoNames.size());

    QStringList list;
    for (unsigned i = 0; i < nGeoObjects; ++i)
        list.append(QString::fromStdString(geoNames[i]));

    if (list.empty())
    {
        this->selectGeoButton->setDisabled(true);
        this->deselectGeoButton->setDisabled(true);
        list.append("[No geometry available.]");
    }
    _allGeo->setStringList(list);
    this->allGeoView->setModel(_allGeo);
    this->selectedGeoView->setModel(_selGeo);
    this->radioAdaptive->toggle(); // default is adaptive meshing
    this->on_radioAdaptive_toggled(true);
}

GMSHPrefsDialog::~GMSHPrefsDialog()
{
    delete _allGeo;
    delete _selGeo;
}

void GMSHPrefsDialog::on_selectGeoButton_pressed()
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

void GMSHPrefsDialog::on_deselectGeoButton_pressed()
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
    if (this->_selGeo->stringList().empty())
    {
        OGSError::box("No geometry selected. Geometric data\n is necessary for mesh generation.");
        return;
    }

    std::vector<std::string> selectedObjects = this->getSelectedObjects(_selGeo->stringList());
    unsigned max_number_of_points_in_quadtree_leaf (10);
    double mesh_density_scaling_pnts(0.5);
    double mesh_density_scaling_stations (0.05);
    double val4(-1);

    if (this->radioAdaptive->isChecked())
    {
        max_number_of_points_in_quadtree_leaf = BaseLib::str2number<unsigned> (
                param1->text().toStdString().c_str());
        if (max_number_of_points_in_quadtree_leaf == 0)
            max_number_of_points_in_quadtree_leaf = 10;
        mesh_density_scaling_pnts = fabs (param2->text().toDouble());
        if (mesh_density_scaling_pnts < sqrt(std::numeric_limits<double>::epsilon()))
            mesh_density_scaling_pnts = 0.5;
        mesh_density_scaling_stations = param3->text().toDouble();
        if (mesh_density_scaling_stations < sqrt(std::numeric_limits<double>::epsilon()))
            mesh_density_scaling_stations = 0.05;
    }
    else
        val4 = param4->text().toDouble();

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
    for (QStringList::iterator it = list.begin(); it != list.end(); ++it)
        indexList.push_back(it->toStdString());
    return indexList;
}
