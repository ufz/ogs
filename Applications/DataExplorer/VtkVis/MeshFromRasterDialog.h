/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-17
 * \brief  Definition of the MeshFromRasterDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ui_MeshFromRaster.h"

#include <string>
#include <QDialog>

namespace MeshLib {
    enum class MeshElemType;
    enum class UseIntensityAs;
}


/**
 * \brief A dialog for specifying the parameters to construct a mesh based on a raster
 */
class MeshFromRasterDialog : public QDialog, private Ui_MeshFromRaster
{
    Q_OBJECT

public:
    /// Constructor
    explicit MeshFromRasterDialog(QDialog* parent = nullptr);
    ~MeshFromRasterDialog() override;

    std::string getMeshName() const { return mesh_name_; }
    std::string getArrayName() const { return array_name_; }
    MeshLib::MeshElemType getElementSelection() const { return element_selection_; }
    MeshLib::UseIntensityAs getIntensitySelection() const { return intensity_selection_; }

private slots:
    void on_elevationButton_toggled(bool isChecked);

    /// Instructions if the OK-Button has been pressed.
    void accept() override;

    /// Instructions if the Cancel-Button has been pressed.
    void reject() override;

private:
    std::string mesh_name_;
    std::string array_name_;
    MeshLib::MeshElemType element_selection_;
    MeshLib::UseIntensityAs intensity_selection_;

};
