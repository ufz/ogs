/**
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (https://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NetCdfConfigureDialog.h"

#include <vtkImageImport.h>

#include <QMessageBox>
#include <QSettings>

#include "GeoLib/Raster.h"
#include "MathLib/Point3d.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/MeshGenerators/RasterToMesh.h"
#include "VtkVis/VtkGeoImageSource.h"
#include "VtkVis/VtkRaster.h"

using namespace netCDF;

// Constructor
NetCdfConfigureDialog::NetCdfConfigureDialog(std::string const& fileName,
                                             QDialog* parent)
    : QDialog(parent),
      _currentFile(fileName.c_str(), NcFile::read),
      _currentMesh(nullptr),
      _currentRaster(nullptr),
      _currentPath(fileName)
{
    setupUi(this);

    int const idx =
        setVariableSelect();  // set up variables of the file in the combobox
    comboBoxVariable->setCurrentIndex(
        idx);  // pre-select the variable with the biggest number of
               // dimensions...valueWithMaxDim()
    _currentVar =
        _currentFile.getVar(comboBoxVariable->itemText(idx).toStdString());

    setDimensionSelect();

    lineEditName->setText(setName());

    this->radioMesh->setChecked(true);
}

NetCdfConfigureDialog::~NetCdfConfigureDialog() = default;

// Instructions if the OK-Button has been pressed.
void NetCdfConfigureDialog::accept()
{
    QMessageBox valueErrorBox;
    if (_currentVar.getDimCount() < 2)
    {
        valueErrorBox.setText("Selected Variable has not enough dimensions.");
        valueErrorBox.exec();
    }
    else if (doubleSpinBoxDim2Start->value() ==
             doubleSpinBoxDim2Start->maximum())
    {
        valueErrorBox.setText("Lon has invalid extend.");
        valueErrorBox.exec();
    }
    else if (doubleSpinBoxDim1Start->value() ==
             doubleSpinBoxDim1Start->maximum())
    {
        valueErrorBox.setText("Lat has invalid extend.");
        valueErrorBox.exec();
    }
    else
    {
        createDataObject();
        this->done(QDialog::Accepted);
    }
}

// Instructions if the Cancel-Button has been pressed.
void NetCdfConfigureDialog::reject()
{
    this->done(QDialog::Rejected);
}

void NetCdfConfigureDialog::on_comboBoxVariable_currentIndexChanged(int /*id*/)
{
    std::string const var_name = comboBoxVariable->currentText().toStdString();
    _currentVar = _currentFile.getVar(var_name);
    setDimensionSelect();
}

// set up x-axis/lat
void NetCdfConfigureDialog::on_comboBoxDim1_currentIndexChanged(int /*id*/)
{
    double firstValue = 0, lastValue = 0;
    unsigned size = 0;
    getDimEdges(comboBoxDim1->currentText().toStdString(), size, firstValue,
                lastValue);
    doubleSpinBoxDim1Start->setValue(firstValue);
    doubleSpinBoxDim1End->setValue(lastValue);
    doubleSpinBoxResolution->setValue(getResolution());
}

// set up y-axis/lon
void NetCdfConfigureDialog::on_comboBoxDim2_currentIndexChanged(int /*id*/)
{
    if (_currentVar.getDimCount() > 1)
    {
        double firstValue = 0, lastValue = 0;
        unsigned size = 0;
        getDimEdges(comboBoxDim2->currentText().toStdString(), size, firstValue,
                    lastValue);
        doubleSpinBoxDim2Start->setValue(firstValue);
        doubleSpinBoxDim2End->setValue(lastValue);
    }
}

// set up time
void NetCdfConfigureDialog::on_comboBoxDim3_currentIndexChanged(int /*id*/)
{
    if (_currentVar.getDimCount() > 2)
    {
        double firstValue = 0, lastValue = 0;
        unsigned size = 0;
        getDimEdges(comboBoxDim3->currentText().toStdString(), size, firstValue,
                    lastValue);
        dateTimeEditDim3->setValue(static_cast<int>(firstValue));
        dateTimeEditDim3->setMinimum(static_cast<int>(firstValue));
        dateTimeEditDim3->setMaximum(static_cast<int>(lastValue));
        lineEditName->setText(setName());
    }
}

// set up additional dimension
void NetCdfConfigureDialog::on_comboBoxDim4_currentIndexChanged(int /*id*/)
{
    if (_currentVar.getDimCount() > 3)
    {
        double firstValue = 0, lastValue = 0;
        unsigned size = 0;
        getDimEdges(comboBoxDim4->currentText().toStdString(), size, firstValue,
                    lastValue);
        spinBoxDim4->setValue(static_cast<int>(firstValue));
        spinBoxDim4->setMinimum(static_cast<int>(firstValue));
        spinBoxDim4->setMaximum(static_cast<int>(lastValue));
    }
}

int NetCdfConfigureDialog::setVariableSelect()
{
    int max_dim = 0;
    int max_dim_idx = 0;
    auto const& names = _currentFile.getVars();
    for (auto [name, var] : names)
    {
        int const var_dim_count = var.getDimCount();
        if (var_dim_count > 1)
        {
            comboBoxVariable->addItem(QString::fromStdString(name));
            if (var_dim_count > max_dim)
            {
                max_dim = var_dim_count;
                max_dim_idx = comboBoxVariable->count() - 1;
            }
        }
    }
    return max_dim_idx;
}

void NetCdfConfigureDialog::setDimensionSelect()
{
    int const dim_count = _currentVar.getDimCount();
    std::array<QComboBox*, 4> dim_box = {
        {comboBoxDim1, comboBoxDim2, comboBoxDim3, comboBoxDim4}};

    for (int i = 0; i < 4; ++i)
    {
        dim_box[i]->clear();
        dim_box[i]->setEnabled(i < dim_count);
    }

    // write dimension-names into selection-boxes
    for (int i = 0; i < dim_count; ++i)
    {
        for (int j = 0; j < dim_count; ++j)
        {
            dim_box[j]->addItem(
                QString::fromStdString(_currentVar.getDim(i).getName()));
        }
    }
    comboBoxDim1->setCurrentIndex(dim_count - 2);
    on_comboBoxDim1_currentIndexChanged(dim_count - 2);
    comboBoxDim2->setCurrentIndex(dim_count - 1);
    on_comboBoxDim2_currentIndexChanged(dim_count - 1);
    // time is only enabled if dim > 2
    dateTimeEditDim3->setEnabled(dim_count > 2);
    // 3rd data dimension is only enabled if dim > 3
    spinBoxDim4->setEnabled(dim_count > 3);

    if (dim_count > 2)
    {
        comboBoxDim3->setCurrentIndex(0);
        on_comboBoxDim3_currentIndexChanged(0);
    }
    else
        dateTimeEditDim3->setSingleStep(0);

    if (dim_count == 4)
    {
        comboBoxDim4->setCurrentIndex(1);
        on_comboBoxDim4_currentIndexChanged(1);
    }
    else
        spinBoxDim4->setValue(0);
}

void NetCdfConfigureDialog::getDimEdges(std::string const& name, unsigned& size,
                                        double& firstValue, double& lastValue)
{
    size = 0;
    firstValue = 0;
    lastValue = 0;
    if (_currentFile.getVar(name).isNull())
        return;

    NcVar const& tmpVarOfDim = _currentFile.getVar(name);
    if ((tmpVarOfDim.getDimCount()) == 1)
    {
        size = tmpVarOfDim.getDim(0).getSize();
        tmpVarOfDim.getVar({0}, {1}, &firstValue);
        tmpVarOfDim.getVar({size - 1}, {1}, &lastValue);
    }
}

int NetCdfConfigureDialog::getTimeStep() const
{
    return dateTimeEditDim3->value();
}

int NetCdfConfigureDialog::getDim4() const
{
    NcVar const& dim3Var =
        _currentFile.getVar(comboBoxDim4->currentText().toStdString());
    std::vector<std::size_t> start{
        static_cast<std::size_t>(spinBoxDim4->value())};
    int value(0);
    dim3Var.getVar(start, {1}, &value);
    if (value < 0)
        value = 0;  // if the value isn't found in the array, set it to 0 as
                    // default...
    return value;
}

double NetCdfConfigureDialog::getResolution()
{
    if (comboBoxDim1->currentIndex() > -1)
    {
        NcVar const& var =
            _currentFile.getVar(comboBoxDim1->currentText().toStdString());
        double firstValue = 0, lastValue = 0;
        unsigned size = 0;
        getDimEdges(var.getName(), size, firstValue, lastValue);
        if (size < 2)
        {
            return 1;
        }

        double interval = fabs(lastValue - firstValue);
        double resolution = (double)interval / (size - 1);
        return resolution;
    }

    return 0;
}

void NetCdfConfigureDialog::createDataObject()
{
    double originLon = 0, originLat = 0;
    double lastLon = 0, lastLat = 0;
    unsigned sizeLon = 0, sizeLat = 0;
    std::string const dim1_name = comboBoxDim1->currentText().toStdString();
    getDimEdges(dim1_name, sizeLat, originLat, lastLat);
    std::string const dim2_name = comboBoxDim2->currentText().toStdString();
    getDimEdges(dim2_name, sizeLon, originLon, lastLon);

    // set up array
    std::vector<double> data_array(sizeLat * sizeLon, 0);

    std::vector<std::size_t> data_origin;
    std::vector<std::size_t> data_length;

    if (_currentVar.getDimCount() > 2)
    {
        // time
        data_origin.push_back(getTimeStep());
        data_length.push_back(1);
        // 3rd dimension
        if (_currentVar.getDimCount() > 3)
        {
            data_origin.push_back(getDim4());
            data_length.push_back(1);
        }
    }

    data_origin.push_back(0);  // x-origin
    data_origin.push_back(0);  // y-origin
    data_length.push_back(sizeLat);
    data_length.push_back(sizeLon);
    _currentVar.getVar(data_origin, data_length, data_array.data());

    std::replace_if(
        data_array.begin(), data_array.end(),
        [](double const& x) { return x <= -9999; }, -9999);

    double origin_x = (originLon < lastLon) ? originLon : lastLon;
    double origin_y = (originLat < lastLat) ? originLat : lastLat;
    MathLib::Point3d origin(std::array<double, 3>{{origin_x, origin_y, 0}});
    double resolution = (doubleSpinBoxResolution->value());

    if (originLat >
        lastLat)  // reverse lines in vertical direction if the original file
                  // has its origin in the northwest corner
        this->reverseNorthSouth(data_array.data(), sizeLon, sizeLat);

    GeoLib::RasterHeader const header = {sizeLon, sizeLat,    1,
                                         origin,  resolution, -9999};
    if (this->radioMesh->isChecked())
    {
        MeshLib::MeshElemType meshElemType = MeshLib::MeshElemType::QUAD;
        MeshLib::UseIntensityAs useIntensity =
            MeshLib::UseIntensityAs::DATAVECTOR;
        if (comboBoxMeshElemType->currentIndex() == 1)
        {
            meshElemType = MeshLib::MeshElemType::TRIANGLE;
        }
        else
        {
            meshElemType = MeshLib::MeshElemType::QUAD;
        }
        if ((comboBoxUseIntensity->currentIndex()) == 1)
        {
            useIntensity = MeshLib::UseIntensityAs::ELEVATION;
        }
        else
        {
            useIntensity = MeshLib::UseIntensityAs::DATAVECTOR;
        }
        _currentMesh = MeshLib::RasterToMesh::convert(
            data_array.data(), header, meshElemType, useIntensity,
            _currentVar.getName());
    }
    else
    {
        vtkImageImport* image =
            VtkRaster::loadImageFromArray(data_array.data(), header);
        _currentRaster = VtkGeoImageSource::New();
        _currentRaster->setImage(image,
                                 QString::fromStdString(this->getName()));
    }
}

QString NetCdfConfigureDialog::setName()
{
    std::string name;
    name.append(_currentPath);
    name.erase(0, name.find_last_of("/") + 1);
    name.erase(name.find_last_of("."));
    return QString::fromStdString(name);
}

std::string NetCdfConfigureDialog::getName()
{
    std::string name = (lineEditName->text()).toStdString();
    QString const date = QString::number(dateTimeEditDim3->value());
    name.append(" - ").append(date.toStdString());
    return name;
}

void NetCdfConfigureDialog::reverseNorthSouth(double* data, std::size_t width,
                                              std::size_t height)
{
    auto* cp_array = new double[width * height];

    for (std::size_t i = 0; i < height; i++)
    {
        for (std::size_t j = 0; j < width; j++)
        {
            std::size_t old_index((width * height) - (width * (i + 1)));
            std::size_t new_index(width * i);
            cp_array[new_index + j] = data[old_index + j];
        }
    }

    std::size_t length(height * width);
    for (std::size_t i = 0; i < length; i++)
        data[i] = cp_array[i];

    delete[] cp_array;
}

void NetCdfConfigureDialog::on_radioMesh_toggled(bool isTrue)
{
    if (isTrue)  // output set to "mesh"
    {
        this->label_2->setEnabled(true);
        this->label_3->setEnabled(false);
        this->comboBoxMeshElemType->setEnabled(true);
        this->comboBoxUseIntensity->setEnabled(false);
    }
    else  // output set to "raster"
    {
        this->label_2->setEnabled(false);
        this->label_3->setEnabled(true);
        this->comboBoxMeshElemType->setEnabled(false);
        this->comboBoxUseIntensity->setEnabled(true);
    }
}
