/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NetCdfConfigureDialog.h"

#include "MathLib/Point3d.h"
#include "GeoLib/Raster.h"
#include "MeshLib/MeshGenerators/RasterToMesh.h"
#include "MeshLib/MeshEnums.h"
#include "VtkVis/VtkGeoImageSource.h"
#include "VtkVis/VtkRaster.h"

#include <QMessageBox>
#include <QSettings>

#include <vtkImageImport.h>

using namespace netCDF;

// Constructor
NetCdfConfigureDialog::NetCdfConfigureDialog(std::string const& fileName, QDialog* parent)
    : QDialog(parent), _currentFile(fileName.c_str(), NcFile::read),
  _currentInitialDateTime(QDateTime()), _currentMesh(nullptr), _currentRaster(nullptr), _currentPath(fileName)
{
    setupUi(this);

    int const idx = setVariableSelect(); // set up variables of the file in the combobox
    comboBoxVariable->setCurrentIndex(idx); //pre-select the variable with the biggest number of dimensions...valueWithMaxDim()
    _currentVar = _currentFile.getVar(comboBoxVariable->itemText(idx).toStdString());

    setDimensionSelect();

    lineEditName->setText(setName());

    this->radioMesh->setChecked(true);
}

NetCdfConfigureDialog::~NetCdfConfigureDialog() = default;

// Instructions if the OK-Button has been pressed.
void NetCdfConfigureDialog::accept()
{
    QMessageBox valueErrorBox;
    if (_currentVar.getDimCount() < 2){
        valueErrorBox.setText("Selected Variable has not enough dimensions.");
        valueErrorBox.exec();
    }else if (doubleSpinBoxDim2Start->value() == doubleSpinBoxDim2Start->maximum()){
        valueErrorBox.setText("Lon has invalid extend.");
        valueErrorBox.exec();
    }else if(doubleSpinBoxDim1Start->value() == doubleSpinBoxDim1Start->maximum()){
        valueErrorBox.setText("Lat has invalid extend.");
        valueErrorBox.exec();
    }else{
        createDataObject();
        this->done(QDialog::Accepted);
    }
}

// Instructions if the Cancel-Button has been pressed.
void NetCdfConfigureDialog::reject()
{
    this->done(QDialog::Rejected);
}

void NetCdfConfigureDialog::on_comboBoxVariable_currentIndexChanged(int id)
{
    std::string const var_name = comboBoxVariable->itemText(_id_map[id]).toStdString();
    _currentVar = _currentFile.getVar(var_name);
    setDimensionSelect();
}

//set up x-axis/lat
void NetCdfConfigureDialog::on_comboBoxDim1_currentIndexChanged(int id)
{
    if (id == -1) id = 0;
    double firstValue=0, lastValue=0;
    unsigned size = 0;
    getDimEdges(id,size,firstValue,lastValue);
    doubleSpinBoxDim1Start->setValue(firstValue);
    doubleSpinBoxDim1End->setValue(lastValue);
    doubleSpinBoxResolution->setValue(getResolution());
}

//set up y-axis/lon
void NetCdfConfigureDialog::on_comboBoxDim2_currentIndexChanged(int id)
{
    if (_currentVar.getDimCount() > 1)
    {
        if (id == -1) id = 0;
        double firstValue=0, lastValue=0;
        unsigned size = 0;
        getDimEdges(id,size,firstValue,lastValue);
        doubleSpinBoxDim2Start->setValue(firstValue);
        doubleSpinBoxDim2End->setValue(lastValue);
    }
}

//set up time
void NetCdfConfigureDialog::on_comboBoxDim3_currentIndexChanged(int id)
{
    if (_currentVar.getDimCount() > 2)
    {
        if (id == -1) id = 0;
        double firstValue=0, lastValue=0;
        unsigned size = 0;
        getDimEdges(id,size,firstValue,lastValue);

        QTime firstTime(0,0,0), lastTime(0,0,0);
        int firstDaysToAdd = 0, lastDaysToAdd = 0;

        getDaysTime(firstValue,firstTime,firstDaysToAdd);
        getDaysTime(lastValue,lastTime,lastDaysToAdd);

        QDate initialDate(1960,1,1);
        QTime initialTime(0,0);

        QDateTime initialDateTime;
        initialDateTime.setDate(initialDate);
        initialDateTime.setTime(initialTime);

        QDateTime firstDateTime = initialDateTime.addDays(firstDaysToAdd);
        firstDateTime.setTime(firstTime);

        QDateTime lastDateTime = initialDateTime.addDays(lastDaysToAdd);
        lastDateTime.setTime(lastTime);

        dateTimeEditDim3->setDateTime(firstDateTime);
        dateTimeEditDim3->setMinimumDateTime(firstDateTime);
        dateTimeEditDim3->setMaximumDateTime(lastDateTime);

        _currentInitialDateTime = initialDateTime;
        lineEditName->setText(setName());
    }
}

//set up additional dimension
void NetCdfConfigureDialog::on_comboBoxDim4_currentIndexChanged(int id)
{
    if (_currentVar.getDimCount() > 3)
    {
        if (id == -1) id = 0;
        double firstValue=0, lastValue=0;
        unsigned size = 0;
        getDimEdges(id,size,firstValue,lastValue);
        // WARNING: Implicit conversion to int in spinBoxDim4->set*()
        spinBoxDim4->setValue(static_cast<int>(firstValue));
        spinBoxDim4->setMinimum(static_cast<int>(firstValue));
        spinBoxDim4->setMaximum(static_cast<int>(lastValue));
    }
}

int NetCdfConfigureDialog::setVariableSelect()
{
    int max_dim = 0;
    int max_dim_idx = 0;
    auto const& names =_currentFile.getVars();
    //auto const n_vars = _currentFile.getVarCount();
    //for (int i = 0; i < n_vars; i++)
    for (auto it = names.cbegin(); it != names.cend(); ++it)
    {
        //NcVar const& focusedVar = _currentFile.getVar(i);
        int const var_dim_count = it->second.getDimCount();
        if (var_dim_count > 1)
        {
            //_id_map.push_back(i);
            comboBoxVariable->addItem(QString::fromStdString(it->first));
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
    int const dim = _currentVar.getDimCount();
    std::array<QComboBox*,4> dim_box = {{ comboBoxDim1, comboBoxDim2, comboBoxDim3, comboBoxDim4 }};

    for (int i = 0; i < 4; ++i)
    {
        dim_box[i]->clear();
        dim_box[i]->setEnabled(i < dim);
    }

    for (int i=0; i < dim; ++i) //write dimension-names into selection-boxes
    {
        for (int j = 0; j < dim; ++j)
            dim_box[j]->addItem(QString::fromStdString(_currentVar.getDim(i).getName()));
    }
    comboBoxDim1->setCurrentIndex(dim-2);
    on_comboBoxDim1_currentIndexChanged(dim-2);
    comboBoxDim2->setCurrentIndex(dim-1);
    on_comboBoxDim2_currentIndexChanged(dim-1);
    dateTimeEditDim3->setEnabled(dim > 2); // time is only enabled if dim > 2
    spinBoxDim4->setEnabled(dim > 3); // add. info is only enabled if dim > 3

    if (dim > 2)
    {
        comboBoxDim3->setCurrentIndex(0);
        on_comboBoxDim3_currentIndexChanged(0);
    }
    else
        dateTimeEditDim3->setDateTime(_currentInitialDateTime);

    if (dim == 4)
    {
        comboBoxDim4->setCurrentIndex(1);
        on_comboBoxDim4_currentIndexChanged(1);
    }
    else
        spinBoxDim4->setValue(0);
}

void NetCdfConfigureDialog::getDimEdges(int dimId, unsigned &size, double &firstValue, double &lastValue)
{
    size = 0;
    firstValue = 0;
    lastValue = 0;
    if (_currentFile.getVar(_currentVar.getDim(dimId).getName()).isNull())
        return;

    NcVar const& tmpVarOfDim = _currentFile.getVar(_currentVar.getDim(dimId).getName());
    if ((tmpVarOfDim.getDimCount()) == 1)
    {
        /*
        int sizeOfDim = tmpVarOfDim.getDim(0).getSize();
        size = sizeOfDim;
        double *arrayOfDimStart = new double[1]; //[1] = {0};
        long edgeOfArray = 1; //[1] = {1};
        long edgeOrigin[1] = {0};
        tmpVarOfDim.set_cur(edgeOrigin);
        tmpVarOfDim->get(arrayOfDimStart,edgeOfArray);
        firstValue = arrayOfDimStart[0];
        double arrayOfDimEnd[1] = {0};
        edgeOrigin[0] = sizeOfDim - 1;
        tmpVarOfDim->set_cur(edgeOrigin);
        tmpVarOfDim->get(arrayOfDimEnd,edgeOfArray);
        lastValue = arrayOfDimEnd[0];
        */
        size = tmpVarOfDim.getDim(0).getSize();
        std::vector<std::size_t> start_val{{0}};
        std::vector<std::size_t> length{{1}};
        tmpVarOfDim.getVar(start_val, length, &firstValue);
        start_val[0] = size-1;
        tmpVarOfDim.getVar(start_val, length, &lastValue);
    }
}

void NetCdfConfigureDialog::getDaysTime(double minSince, QTime &time, int &days)
{
    auto tmpMin = (long)minSince;
    long minuits = tmpMin % 60;
    long tmpHours = tmpMin / 60;
    long hours = tmpHours % 24;
    days = (int)(tmpHours / 24);
    time.setHMS(hours, minuits, 0);
}

long NetCdfConfigureDialog::convertDateToMinutes(QDateTime initialDateTime, QDate selectedDate, QTime selectedTime)
{
    long tmpInitialToSelectedDate = static_cast<long>(selectedDate.daysTo(initialDateTime.date()));
    long selectedDays = - tmpInitialToSelectedDate * 24 * 60;
    long selectedMinutes = (selectedTime.hour() * 60) + selectedTime.minute() + selectedDays;
    return selectedMinutes;
}

int NetCdfConfigureDialog::getTimeStep()
{
    NcVar const& timeVar =
        _currentFile.getVar(comboBoxDim2->currentText().toStdString());

    double const datesToMinutes = convertDateToMinutes(_currentInitialDateTime,dateTimeEditDim3->date(),dateTimeEditDim3->time());

    double timeArray[1] = {datesToMinutes};
    /*
    double currentTime = timeVar.get_index(timeArray);
    if (currentTime < 0) currentTime=0; //if the value isn't found in the array, set it to 0 as default...
    return static_cast<int>(currentTime);
    */
    /** TODO **/
    return 0;
}

int NetCdfConfigureDialog::getDim4() const
{
    /** TODO **/
    /*
    NcVar const& dim3Var =
        _currentFile.getVar(comboBoxDim4->currentText().toStdString());
    int timeArray[1] = {spinBoxDim4->value()};
    int currentValueDim3 = dim3Var.get_index(timeArray);
    if (currentValueDim3 < 0)
        currentValueDim3=0; //if the value isn't found in the array, set it to 0 as default...
    return currentValueDim3;
    */
    return 0;
}

double NetCdfConfigureDialog::getResolution()
{
    if (comboBoxDim1->currentIndex() > -1)
    {
        NcVar const& latVar =
            _currentFile.getVar(comboBoxDim1->currentText().toStdString());
        double firstValue = 0, lastValue = 0;
        unsigned size = 0;
        getDimEdges(latVar.getId(), size, firstValue, lastValue);
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
    auto* length = new std::size_t[_currentVar.getDimCount()];
    double originLon = 0, originLat = 0;
    double lastLon = 0, lastLat = 0;
    unsigned sizeLon = 0, sizeLat = 0;
    getDimEdges(comboBoxDim1->currentIndex(), sizeLat, originLat, lastLat);
    getDimEdges(comboBoxDim2->currentIndex(), sizeLon, originLon, lastLon);


    for(int i=0; i < _currentVar.getDimCount(); i++) length[i]=1;

    // set array edges: lat x lon
    length[comboBoxDim1->currentIndex()]=sizeLat;
    length[comboBoxDim2->currentIndex()]=sizeLon;

    // set up array
    auto* data_array = new double[sizeLat * sizeLon];
    for(std::size_t i=0; i < (sizeLat*sizeLon); i++) data_array[i]=0;

    //Time-Dimension:
    /*
    if (_currentVar.getDimCount() > 2)
    {
        auto* newOrigin = new long[_currentVar.getDimCount()];
        for (int i=0; i < _currentVar.getDimCount(); i++) newOrigin[i]=0;
        newOrigin[comboBoxDim3->currentIndex()] = getTimeStep(); //set origin to selected time
        _currentVar.set_cur(newOrigin);
        //Dimension 4:
        if (_currentVar.getDimCount() > 3) newOrigin[comboBoxDim4->currentIndex()] = getDim4(); //if there are is a 4th dimension
        delete [] newOrigin;
    }
    */
    std::vector<std::size_t> data_origin {{0, 0}};
    std::vector<std::size_t> data_length {{sizeLat, sizeLon}};
    _currentVar.getVar(data_origin, data_length, data_array);  // create Array of Values

    for (std::size_t i=0; i < (sizeLat*sizeLon); i++)
    {
        //data_array[i] = data_array[i] - 273; // convert from kalvin to celsius
        if (data_array[i] < -9999 ) data_array[i] = -9999; // all values < -10000, set to "no-value"
    }

    double origin_x = (originLon < lastLon) ? originLon : lastLon;
    double origin_y = (originLat < lastLat) ? originLat : lastLat;
    MathLib::Point3d origin(std::array<double,3>{{origin_x, origin_y, 0}});
    double resolution = (doubleSpinBoxResolution->value());

    if (originLat > lastLat) // reverse lines in vertical direction if the original file has its origin in the northwest corner
        this->reverseNorthSouth(data_array, sizeLon, sizeLat);

    GeoLib::RasterHeader const header = {sizeLon, sizeLat,    1,
                                         origin,  resolution, -9999};
    if (this->radioMesh->isChecked())
    {
        MeshLib::MeshElemType meshElemType = MeshLib::MeshElemType::QUAD;
        MeshLib::UseIntensityAs useIntensity = MeshLib::UseIntensityAs::DATAVECTOR;
        if (comboBoxMeshElemType->currentIndex() == 1)
        {
            meshElemType = MeshLib::MeshElemType::TRIANGLE;
        }else{
            meshElemType = MeshLib::MeshElemType::QUAD;
        }
        if ((comboBoxUseIntensity->currentIndex()) == 1)
        {
            useIntensity = MeshLib::UseIntensityAs::ELEVATION;
        }else{
            useIntensity = MeshLib::UseIntensityAs::DATAVECTOR;
        }
        _currentMesh = MeshLib::RasterToMesh::convert(
            data_array, header, meshElemType, useIntensity, _currentVar.getName());
    }
    else
    {
        vtkImageImport* image = VtkRaster::loadImageFromArray(data_array, header);
        _currentRaster = VtkGeoImageSource::New();
        _currentRaster->setImage(image, QString::fromStdString(this->getName()), origin[0], origin[1], resolution);
    }

    delete[] length;
    delete[] data_array;
}

QString NetCdfConfigureDialog::setName()
{
    std::string name;
    name.append(_currentPath);
    name.erase(0,name.find_last_of("/")+1);
    name.erase(name.find_last_of("."));
    return QString::fromStdString(name);
}

std::string NetCdfConfigureDialog::getName()
{
    std::string name = (lineEditName->text()).toStdString();
    QString date = dateTimeEditDim3->date().toString(Qt::LocalDate);
    name.append(" - ").append(date.toStdString());
    return name;
}

void NetCdfConfigureDialog::reverseNorthSouth(double* data, std::size_t width, std::size_t height)
{
    auto* cp_array = new double[width * height];

    for (std::size_t i=0; i<height; i++)
    {
        for (std::size_t j=0; j<width; j++)
        {
            std::size_t old_index((width*height)-(width*(i+1)));
            std::size_t new_index(width*i);
            cp_array[new_index+j] = data[old_index+j];
        }
    }

    std::size_t length(height*width);
    for (std::size_t i=0; i<length; i++)
        data[i] = cp_array[i];

    delete[] cp_array;
}

void NetCdfConfigureDialog::on_radioMesh_toggled(bool isTrue)
{
    if (isTrue) // output set to "mesh"
    {
        this->label_2->setEnabled(true);
        this->label_3->setEnabled(false);
        this->comboBoxMeshElemType->setEnabled(true);
        this->comboBoxUseIntensity->setEnabled(false);
    }
    else // output set to "raster"
    {
        this->label_2->setEnabled(false);
        this->label_3->setEnabled(true);
        this->comboBoxMeshElemType->setEnabled(false);
        this->comboBoxUseIntensity->setEnabled(true);
    }
}
