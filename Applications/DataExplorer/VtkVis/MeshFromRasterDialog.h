// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
