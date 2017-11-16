/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-17
 * \brief  Definition of the MeshFromRasterDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
    MeshFromRasterDialog(QDialog* parent = nullptr);
    ~MeshFromRasterDialog(void) override;

    std::string getMeshName() const { return _mesh_name; }
    std::string getArrayName() const { return _array_name; }
    MeshLib::MeshElemType getElementSelection() const { return _element_selection; }
    MeshLib::UseIntensityAs getIntensitySelection() const { return _intensity_selection; }

private slots:
    void on_elevationButton_toggled(bool isChecked);

    /// Instructions if the OK-Button has been pressed.
    void accept() override;

    /// Instructions if the Cancel-Button has been pressed.
    void reject() override;

private:
    std::string _mesh_name;
    std::string _array_name;
    MeshLib::MeshElemType _element_selection;
    MeshLib::UseIntensityAs _intensity_selection;

};
