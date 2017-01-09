/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NETCDFCONFIGUREDIALOG_H
#define NETCDFCONFIGUREDIALOG_H

#include <netcdfcpp.h>

#include <QDialog>
#include "ui_NetCdfConfigure.h"

namespace MeshLib {
    class Mesh;
}

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
    NetCdfConfigureDialog(const std::string &fileName, QDialog* parent = 0);
    ~NetCdfConfigureDialog(void);
    MeshLib::Mesh* getMesh() { return _currentMesh; };
    std::string getName();
    VtkGeoImageSource* getRaster() { return _currentRaster; };

private slots:
    void accept();
    void reject();
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
    void getDimEdges(int dimId, unsigned &size, double &firstValue, double &lastValue);
    void getDaysTime(double minSince, QTime &time, int &days);
    long convertDateToMinutes(QDateTime initialDateTime,QDate selectedDate, QTime selectedTime);
    void createDataObject();
    int getTimeStep();
    int getDim4() const;
    double getResolution();
    QString setName();
    void reverseNorthSouth(double* data, std::size_t width, std::size_t height);

    NcFile *_currentFile;
    NcVar *_currentVar;
    QDateTime _currentInitialDateTime;
    MeshLib::Mesh* _currentMesh;
    VtkGeoImageSource* _currentRaster;
    std::string _currentPath;
    std::vector<int> _id_map;
};

#endif //NETCDFCONFIGUREDIALOG_H

