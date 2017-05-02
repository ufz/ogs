/**
 * \file
 * \author Karsten Rink
 * \date   2013-05-29
 * \brief  Definition of the MergeGeometriesDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ui_MergeGeometries.h"
#include <QDialog>

class QStringListModel;

namespace GeoLib
{
class GEOObjects;
}

/**
 * \brief A dialog window for setting preferences for GMSH
 */
class MergeGeometriesDialog : public QDialog, private Ui_MergeGeometries
{
    Q_OBJECT

public:
    MergeGeometriesDialog(GeoLib::GEOObjects& geoObjects,
                          QDialog* parent = nullptr);
    ~MergeGeometriesDialog(void) override;

    /// Returns a vector of selected geometries
    std::vector<std::string> const getSelectedGeometries() const;

    /// Returns the name of the new merged geometry
    std::string getGeometryName() const;

private:
    GeoLib::GEOObjects& _geo_objects;
    QStringListModel* _allGeo;
    QStringListModel* _selGeo;

private slots:
    void on_selectGeoButton_pressed();
    void on_deselectGeoButton_pressed();

    /// Instructions if the OK-Button has been pressed.
    void accept() override;

    /// Instructions if the Cancel-Button has been pressed.
    void reject() override;
};
