//file NetCDFConfigureDialog.cpp
//CH Initial implementation

#include "NetCdfConfigureDialog.h"

#include "VtkMeshConverter.h"
#include "GridAdapter.h"
#include "VtkGeoImageSource.h"
#include "VtkRaster.h"

#include <QMessageBox>
#include <QSettings>

#include <vtkImageImport.h>

// Constructor
NetCdfConfigureDialog::NetCdfConfigureDialog(const std::string &fileName, QDialog* parent) 
	: QDialog(parent), _currentFile(new NcFile(fileName.c_str(), NcFile::ReadOnly)), 
	  _currentInitialDateTime(QDateTime()), _currentMesh(NULL), _currentRaster(NULL), _currentPath(fileName)
{
	setupUi(this);

	setVariableSelect(); // set up variables of the file in the combobox
	comboBoxVariable->setCurrentIndex(valueWithMaxDim()); //pre-select the variable with the biggest number of dimensions...valueWithMaxDim()
	
	_currentVar = _currentFile->get_var(comboBoxVariable->currentIndex());

	setDimensionSelect();

	lineEditName->setText(setName());

	this->radioMesh->setChecked(true);

}

NetCdfConfigureDialog::~NetCdfConfigureDialog()
{
}

// Instructions if the OK-Button has been pressed.
void NetCdfConfigureDialog::accept()
{
	QMessageBox valueErrorBox;
	if (_currentVar->num_dims() < 3){
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
		delete _currentFile;
		this->done(QDialog::Accepted);
	}
}

// Instructions if the Cancel-Button has been pressed.
void NetCdfConfigureDialog::reject()
{
	delete _currentFile;
	this->done(QDialog::Rejected);
}

int NetCdfConfigureDialog::valueWithMaxDim()
{
	int idMaxDim = 0;
	for (int i=0; i < _currentFile->num_dims(); i++)
	{
		if ((_currentFile->get_var(i)->num_dims()) > 0) idMaxDim = i + 1;
	}
	return idMaxDim;
}

void NetCdfConfigureDialog::on_comboBoxVariable_currentIndexChanged(int id)
{
	_currentVar = _currentFile->get_var(id);
	setDimensionSelect();
}

//set up x-axis/lat
void NetCdfConfigureDialog::on_comboBoxDim1_currentIndexChanged(int id) 
{
	if (id == -1) id = 0;
	double firstValue=0, lastValue=0;
	size_t size = 0;
	getDimEdges(id,size,firstValue,lastValue);
	doubleSpinBoxDim1Start->setValue(firstValue);
	doubleSpinBoxDim1End->setValue(lastValue);
	doubleSpinBoxResolution->setValue(getResolution());
}

//set up y-axis/lon
void NetCdfConfigureDialog::on_comboBoxDim2_currentIndexChanged(int id)
{
	if (_currentVar->num_dims() > 1)
	{
		if (id == -1) id = 0;
		double firstValue=0, lastValue=0;
		size_t size = 0;
		getDimEdges(id,size,firstValue,lastValue);
		doubleSpinBoxDim2Start->setValue(firstValue);
		doubleSpinBoxDim2End->setValue(lastValue);
	}
}

//set up time
void NetCdfConfigureDialog::on_comboBoxDim3_currentIndexChanged(int id)
{
	if (_currentVar->num_dims() > 2)
	{
		if (id == -1) id = 0;
		double firstValue=0, lastValue=0;
		size_t size = 0;
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
	if (_currentVar->num_dims() > 3)
	{
		if (id == -1) id = 0;
		double firstValue=0, lastValue=0;
		size_t size = 0;
		getDimEdges(id,size,firstValue,lastValue);
		// WARNING: Implicit conversion to int in spinBoxDim4->set*()
		spinBoxDim4->setValue(firstValue);
		spinBoxDim4->setMinimum(firstValue);
		spinBoxDim4->setMaximum(lastValue);
	}
}

void NetCdfConfigureDialog::setVariableSelect()
{
	for (int i=0; i<(_currentFile->num_vars()); i++)
	{
		NcVar *focusedVar = _currentFile->get_var(i);
		if (focusedVar->num_dims() > 0) comboBoxVariable->addItem(focusedVar->name());
	}
}

void NetCdfConfigureDialog::setDimensionSelect()
{
	comboBoxDim1->clear();
	comboBoxDim2->clear();
	comboBoxDim3->clear();
	comboBoxDim4->clear();
	for (int i=0; i < _currentVar->num_dims(); i++) //write dimension-names into selection-boxes
	{
		comboBoxDim1->addItem(_currentVar->get_dim(i)->name());
		comboBoxDim2->addItem(_currentVar->get_dim(i)->name());
		comboBoxDim3->addItem(_currentVar->get_dim(i)->name());
		comboBoxDim4->addItem(_currentVar->get_dim(i)->name());
	}
	if (_currentVar->num_dims() < 4) 
	{
		comboBoxDim4->setEnabled(false);comboBoxDim4->clear();
		spinBoxDim4->setEnabled(false);spinBoxDim4->setValue(0);
	}else{  //pre-set dimension selection, typical for 4-d-nc-files:
		comboBoxDim4->setEnabled(true);
		spinBoxDim4->setEnabled(true);
		comboBoxDim1->setCurrentIndex(2);on_comboBoxDim1_currentIndexChanged(2);
		comboBoxDim2->setCurrentIndex(3);on_comboBoxDim2_currentIndexChanged(3);
		comboBoxDim3->setCurrentIndex(0);on_comboBoxDim3_currentIndexChanged(0);
		comboBoxDim4->setCurrentIndex(1);on_comboBoxDim4_currentIndexChanged(1);
	}
	if (_currentVar->num_dims() == 3) //pre-set dimension selection, typical for 3-d-nc-files:
	{
		comboBoxDim1->setCurrentIndex(1);on_comboBoxDim1_currentIndexChanged(1);
		comboBoxDim2->setCurrentIndex(2);on_comboBoxDim2_currentIndexChanged(2);
		comboBoxDim3->setCurrentIndex(0);on_comboBoxDim3_currentIndexChanged(0);
	}
	if (_currentVar->num_dims() < 3)
	{
		comboBoxDim3->setEnabled(false);comboBoxDim3->clear();
		dateTimeEditDim3->setEnabled(false);dateTimeEditDim3->setDateTime(_currentInitialDateTime);

	}else{
		comboBoxDim3->setEnabled(true);
		dateTimeEditDim3->setEnabled(true);
	}
	if (_currentVar->num_dims() < 2) 
	{
		comboBoxDim2->setEnabled(false);comboBoxDim2->clear();
		doubleSpinBoxDim2Start->setValue(0); doubleSpinBoxDim2End->setValue(0);
	}else{
		comboBoxDim2->setEnabled(true);
	}

}

void NetCdfConfigureDialog::getDimEdges(int dimId, size_t &size, double &firstValue, double &lastValue)
{
	if ((_currentFile->get_var(_currentVar->get_dim(dimId)->name())) != NULL) 
	{
		NcVar *tmpVarOfDim = _currentFile->get_var(_currentVar->get_dim(dimId)->name());
		if ((tmpVarOfDim->num_dims()) == 1)
		{
			int sizeOfDim = tmpVarOfDim->get_dim(0)->size();
			size = sizeOfDim;
			double arrayOfDimStart[1] = {0};
#ifdef VTK_NETCDF_FOUND
			size_t edgeOfArray[1] = {1};
#else
			long edgeOfArray[1] = {1};
#endif
			long edgeOrigin[1] = {0};
			tmpVarOfDim->set_cur(edgeOrigin);
			tmpVarOfDim->get(arrayOfDimStart,edgeOfArray);
			firstValue = arrayOfDimStart[0];
			double arrayOfDimEnd[1] = {0};
			edgeOrigin[0] = sizeOfDim - 1;
			tmpVarOfDim->set_cur(edgeOrigin);
			tmpVarOfDim->get(arrayOfDimEnd,edgeOfArray);
			lastValue = arrayOfDimEnd[0];
		}
	}else{
		size = 0;
		firstValue = 0;
		lastValue = 0;
	}
}

void NetCdfConfigureDialog::getDaysTime(double minSince, QTime &time, int &days)
{
	long tmpMin = (long) minSince;
	long minuits = tmpMin % 60;
	long tmpHours = tmpMin / 60;
	long hours = tmpHours % 24;
	days = (int) (tmpHours / 24);
	time.setHMS(hours,minuits,0);
}

long NetCdfConfigureDialog::convertDateToMinutes(QDateTime initialDateTime, QDate selectedDate, QTime selectedTime)
{
	int tmpInitialToSelectedDate = (selectedDate.daysTo(initialDateTime.date()));
	long selectedDays = - tmpInitialToSelectedDate * 24 * 60;
	long selectedMinutes = (selectedTime.hour() * 60) + selectedTime.minute() + selectedDays;
	return selectedMinutes;
}

int NetCdfConfigureDialog::getTimeStep()
{
	NcVar* timeVar = _currentFile->get_var(comboBoxDim2->currentIndex());
	
	const double datesToMinutes = convertDateToMinutes(_currentInitialDateTime,dateTimeEditDim3->date(),dateTimeEditDim3->time());

	double timeArray[1] = {datesToMinutes};
	double currentTime = timeVar->get_index(timeArray);
	if (currentTime < 0) currentTime=0; //if the value isn't found in the array, set it to 0 as default...
	return currentTime;
}	

int NetCdfConfigureDialog::getDim4()
{
	NcVar* dim3Var = _currentFile->get_var(comboBoxDim4->currentIndex());
	double timeArray[1] = {static_cast<double>(spinBoxDim4->value())};
	double currentValueDim3 = dim3Var->get_index(timeArray);
	if (currentValueDim3 < 0) currentValueDim3=0; //if the value isn't found in the array, set it to 0 as default...
	return currentValueDim3;
}	

double NetCdfConfigureDialog::getResolution()
{
	if (comboBoxDim1->currentIndex() > -1)
	{
		NcVar *latVar = _currentFile->get_var(comboBoxDim1->currentIndex());
		double firstValue=0, lastValue=0;
		size_t size=0;
		getDimEdges(latVar->id(),size,firstValue,lastValue);
		if (size < 2)
		{
			return 1;
		}
		else
		{
			double interval = fabs(lastValue-firstValue);
			double resolution = (double)interval/(size-1);
			return resolution;
		}
	}
	else
	{
		return 0;
	}
}

void NetCdfConfigureDialog::createDataObject()
{
#ifdef VTK_NETCDF_FOUND
	size_t* length = new size_t[_currentVar->num_dims()];
#else
	long* length = new long[_currentVar->num_dims()];
#endif
	double originLon = 0, originLat = 0;
	double lastLon = 0, lastLat = 0;
	size_t sizeLon = 0, sizeLat = 0;
	getDimEdges(comboBoxDim1->currentIndex(), sizeLat, originLat, lastLat);
	getDimEdges(comboBoxDim2->currentIndex(), sizeLon, originLon, lastLon);

	for(int i=0; i < _currentVar->num_dims(); i++) length[i]=1;

	// set array edges: lat x lon
	length[comboBoxDim1->currentIndex()]=sizeLat;
	length[comboBoxDim2->currentIndex()]=sizeLon;

	// set up array
	double* data_array = new double[sizeLat*sizeLon];
	for(size_t i=0; i < (sizeLat*sizeLon); i++) data_array[i]=0;

	//Time-Dimension:
	if (_currentVar->num_dims() > 2)
	{
		long* newOrigin = new long[_currentVar->num_dims()];
		for (int i=0; i < _currentVar->num_dims(); i++) newOrigin[i]=0;
		newOrigin[comboBoxDim3->currentIndex()] = getTimeStep(); //set origin to selected time
		_currentVar->set_cur(newOrigin);
		//Dimension 4:
		if (_currentVar->num_dims() > 3) newOrigin[comboBoxDim4->currentIndex()] = getDim4(); //if there are is a 4th dimension 
		delete newOrigin;
	}
	
	_currentVar->get(data_array,length); //create Array of Values

	for (size_t i=0; i < (sizeLat*sizeLon); i++)
	{
		//data_array[i] = data_array[i] - 273; // convert from kalvin to celsius
		if (data_array[i] < -9999 ) data_array[i] = -9999; // all values < -10000, set to "no-value"
	}
		
	double origin_x = (originLon < lastLon) ? originLon : lastLon;
	double origin_y = (originLat < lastLat) ? originLat : lastLat;
	double originNetCdf[3] = {origin_x, origin_y, 0};

	double resolution = (doubleSpinBoxResolution->value());

	if (originLat > lastLat) // reverse lines in vertical direction if the original file has its origin in the northwest corner
		this->reverseNorthSouth(data_array, sizeLon, sizeLat);

	if (this->radioMesh->isChecked())
	{
		MshElemType::type meshElemType = MshElemType::QUAD;
		UseIntensityAs::type useIntensity = UseIntensityAs::MATERIAL;
		if ((comboBoxMeshElemType->currentIndex()) == 1) 
		{
			meshElemType = MshElemType::TRIANGLE;
		}else{
			meshElemType = MshElemType::QUAD;
		}
		if ((comboBoxUseIntensity->currentIndex()) == 1)
		{
			useIntensity = UseIntensityAs::ELEVATION;
		}else{
			useIntensity = UseIntensityAs::MATERIAL;
		}

		_currentMesh = VtkMeshConverter::convertImgToMesh(data_array,originNetCdf,sizeLon,sizeLat,resolution,meshElemType,useIntensity);
	}
	else
	{
		vtkImageImport* image = VtkRaster::loadImageFromArray(data_array, originNetCdf[0], originNetCdf[1], sizeLon, sizeLat, resolution, -9999.0);
		_currentRaster = VtkGeoImageSource::New();
		_currentRaster->setImage(image, QString::fromStdString(this->getName()), originNetCdf[0], originNetCdf[1], resolution);
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
	QString qstr = QString::fromStdString(name);
	return qstr;
}

std::string NetCdfConfigureDialog::getName()
{
	std::string name = (lineEditName->text()).toStdString();
	QString date = dateTimeEditDim3->date().toString(Qt::LocalDate);
	name.append(" - ").append(date.toStdString());
	return name;	
}

void NetCdfConfigureDialog::reverseNorthSouth(double* data, size_t width, size_t height)
{
	double* cp_array = new double[width*height];

	for (size_t i=0; i<height; i++)
	{
		for (size_t j=0; j<width; j++)
		{
			size_t old_index((width*height)-(width*(i+1)));
			size_t new_index(width*i);
			cp_array[new_index+j] = data[old_index+j];
		}
	}

	size_t length(height*width);
	for (size_t i=0; i<length; i++)
		data[i] = cp_array[i];

	delete[] cp_array;
}

void NetCdfConfigureDialog::on_radioMesh_toggled(bool isTrue)
{
	if (isTrue) // output set to "mesh"
	{
		this->label_2->setEnabled(true);
		this->label_3->setEnabled(true);
		this->comboBoxMeshElemType->setEnabled(true);
		this->comboBoxUseIntensity->setEnabled(true);
	}
	else // output set to "raster"
	{
		this->label_2->setEnabled(false);
		this->label_3->setEnabled(false);
		this->comboBoxMeshElemType->setEnabled(false);
		this->comboBoxUseIntensity->setEnabled(false);
	}
}


