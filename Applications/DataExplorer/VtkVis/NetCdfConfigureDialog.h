/**
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (https://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <QDialog>
#include <netcdf>

#include "MeshLib/Mesh.h"
#include "ui_NetCdfConfigure.h"

class VtkGeoImageSource;

/**
 * \brief The dialog for converting data from NetCDF-files into OGS data structures.
 * While NetCDF files can include data ranging from dimensionality 0 (scalars) to 4
 * (e.g. 3d arrays over time), only arrays of dimensionality 2 or higher can be
 * imported as this is the minimum requirement for a meaningful conversion into
 * raster- or mesh-data-objects. Scalars or vector-variables are not selectable from
 * the selection menu.
 */
class NetCdfConfigureDialog : public QDialog, private Ui_NetCdfConfigure
{
    Q_OBJECT

public:
    explicit NetCdfConfigureDialog(const std::string& fileName,
                                   QDialog* parent = nullptr);
    ~NetCdfConfigureDialog() override;
    MeshLib::Mesh* getMesh() { return _currentMesh.get(); };
    std::string getName();
    VtkGeoImageSource* getRaster() { return _currentRaster; };

private slots:
    void accept() override;
    void reject() override;
    void on_comboBoxVariable_currentIndexChanged(int id);
    void on_comboBoxDim1_currentIndexChanged(int id);
    void on_comboBoxDim2_currentIndexChanged(int id);
    void on_comboBoxDim3_currentIndexChanged(int id);
    void on_comboBoxDim4_currentIndexChanged(int id);
    void on_radioMesh_toggled(bool isTrue);

private:
    /// Fills the combobox with all applicable variables and
    /// returns the index of the first variable with the highest dimension.
    int setVariableSelect();
    void setDimensionSelect();
    void getDimEdges(std::string const& name,
                     unsigned& size,
                     double& firstValue,
                     double& lastValue);
    void createDataObject();
    int getTimeStep() const;
    int getDim4() const;
    double getResolution();
    QString setName();
    void reverseNorthSouth(double* data, std::size_t width, std::size_t height);

    netCDF::NcFile _currentFile;
    netCDF::NcVar _currentVar;
    std::unique_ptr<MeshLib::Mesh> _currentMesh;
    VtkGeoImageSource* _currentRaster;
    std::string _currentPath;
};
