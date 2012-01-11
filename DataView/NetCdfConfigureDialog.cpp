//file NetCDFConfigureDialog.cpp
//CH Initial implementation

#include "NetCdfConfigureDialog.h"

// Constructor
NetCdfConfigureDialog::NetCdfConfigureDialog(char* fileName, QDialog* parent) : QDialog(parent)
{
	setupUi(this);
	currentFile = new NcFile(fileName,NcFile::ReadOnly);
	currentPath = fileName;

	setVariableSelect(); // set up variables of the file in the combobox
	comboBoxVariable->setCurrentIndex(valueWithMaxDim()); //pre-select the variable with the biggest number of dimensions...valueWithMaxDim()
	
	currentVar = currentFile->get_var(comboBoxVariable->currentIndex());

	setDimensionSelect();

	lineEditName->setText(setName());

}

NetCdfConfigureDialog::~NetCdfConfigureDialog()
{
}

// Instructions if the OK-Button has been pressed.
void NetCdfConfigureDialog::accept()
{
	QMessageBox valueErrorBox;
	if (currentVar->num_dims() < 3){
		valueErrorBox.setText("Selected Variable has less dimensions.");
		valueErrorBox.exec();
	}else if (doubleSpinBoxDim2Start->value() == doubleSpinBoxDim2Start->maximum()){
		valueErrorBox.setText("Lon has invalid dimension.");
		valueErrorBox.exec();
	}else if(doubleSpinBoxDim1Start->value() == doubleSpinBoxDim1Start->maximum()){
		valueErrorBox.setText("Lat has invalid dimension.");
		valueErrorBox.exec();
	}else{
		createMesh();
		delete currentFile;
		this->done(QDialog::Accepted);
	}
}

// Instructions if the Cancel-Button has been pressed.
void NetCdfConfigureDialog::reject()
{
	currentMesh = NULL;
	delete currentFile;
	this->done(QDialog::Rejected);
}

int NetCdfConfigureDialog::valueWithMaxDim()
{
	int idMaxDim = 0;
	for (int i=0; i < currentFile->num_dims(); i++)
	{
		if ((currentFile->get_var(i)->num_dims()) > 0) idMaxDim = i + 1;
	}
	return idMaxDim;
}

void NetCdfConfigureDialog::on_comboBoxVariable_currentIndexChanged(int id)
{
	currentVar = currentFile->get_var(id);
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
	if (currentVar->num_dims() > 1)
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
	if (currentVar->num_dims() > 2)
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

		currentInitialDateTime = initialDateTime;
		lineEditName->setText(setName());
	}
}

//set up additional dimension
void NetCdfConfigureDialog::on_comboBoxDim4_currentIndexChanged(int id)
{
	if (currentVar->num_dims() > 3)
	{
		if (id == -1) id = 0;
		double firstValue=0, lastValue=0;
		size_t size = 0;
		getDimEdges(id,size,firstValue,lastValue);
		spinBoxDim4->setValue(firstValue);
		spinBoxDim4->setMinimum(firstValue);
		spinBoxDim4->setMaximum(lastValue);
	}
}

void NetCdfConfigureDialog::setVariableSelect()
{
	for (int i=0; i<(currentFile->num_vars()); i++)
	{
		NcVar *focusedVar = currentFile->get_var(i);
		if (focusedVar->num_dims() > 0) comboBoxVariable->addItem(focusedVar->name());
	}
}

void NetCdfConfigureDialog::setDimensionSelect()
{
	comboBoxDim1->clear();
	comboBoxDim2->clear();
	comboBoxDim3->clear();
	comboBoxDim4->clear();
	for (int i=0; i < currentVar->num_dims(); i++) //write dimension-names into selection-boxes
	{
		comboBoxDim1->addItem(currentVar->get_dim(i)->name());
		comboBoxDim2->addItem(currentVar->get_dim(i)->name());
		comboBoxDim3->addItem(currentVar->get_dim(i)->name());
		comboBoxDim4->addItem(currentVar->get_dim(i)->name());
	}
	if (currentVar->num_dims() < 4) 
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
	if (currentVar->num_dims() == 3) //pre-set dimension selection, typical for 3-d-nc-files:
	{
		comboBoxDim1->setCurrentIndex(1);on_comboBoxDim1_currentIndexChanged(1);
		comboBoxDim2->setCurrentIndex(2);on_comboBoxDim2_currentIndexChanged(2);
		comboBoxDim3->setCurrentIndex(0);on_comboBoxDim3_currentIndexChanged(0);
	}
	if (currentVar->num_dims() < 3)
	{
		comboBoxDim3->setEnabled(false);comboBoxDim3->clear();
		dateTimeEditDim3->setEnabled(false);dateTimeEditDim3->setDateTime(currentInitialDateTime);

	}else{
		comboBoxDim3->setEnabled(true);
		dateTimeEditDim3->setEnabled(true);
	}
	if (currentVar->num_dims() < 2) 
	{
		comboBoxDim2->setEnabled(false);comboBoxDim2->clear();
		doubleSpinBoxDim2Start->setValue(0); doubleSpinBoxDim2End->setValue(0);
	}else{
		comboBoxDim2->setEnabled(true);
	}

}

void NetCdfConfigureDialog::getDimEdges(int dimId, size_t &size, double &firstValue, double &lastValue)
{
	if ((currentFile->get_var(currentVar->get_dim(dimId)->name())) != NULL) 
	{
		NcVar *tmpVarOfDim = currentFile->get_var(currentVar->get_dim(dimId)->name());
		if ((tmpVarOfDim->num_dims()) == 1)
		{
			int sizeOfDim = tmpVarOfDim->get_dim(0)->size();
			size = sizeOfDim;
			double arrayOfDimStart[1] = {0};
			size_t edgeOfArray[1] = {1};
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
	QTime selectedTime = dateTimeEditDim3->time();
	QDate selectedDate = dateTimeEditDim3->date();

	NcVar* timeVar = currentFile->get_var(comboBoxDim2->currentIndex());
	
	int datesToMinutes = convertDateToMinutes(currentInitialDateTime,dateTimeEditDim3->date(),dateTimeEditDim3->time());

	double timeArray[1] = {datesToMinutes};
	double currentTime = timeVar->get_index(timeArray);
	if (currentTime < 0) currentTime=0; //if the value isn't found in the array, set it to 0 as default...
	return currentTime;
}	

int NetCdfConfigureDialog::getDim4()
{
	NcVar* dim3Var = currentFile->get_var(comboBoxDim4->currentIndex());
	double timeArray[1] = {spinBoxDim4->value()};
	double currentValueDim3 = dim3Var->get_index(timeArray);
	if (currentValueDim3 < 0) currentValueDim3=0; //if the value isn't found in the array, set it to 0 as default...
	return currentValueDim3;
}	

int NetCdfConfigureDialog::getResolution()
{
	if (comboBoxDim1->currentIndex() > -1)
	{
		NcVar *latVar = currentFile->get_var(comboBoxDim1->currentIndex());
		double firstValue=0, lastValue=0;
		size_t size=0;
		getDimEdges(latVar->id(),size,firstValue,lastValue);
		if (size < 2)
		{
			return 1;
		}else{
			int resolution = ((lastValue-firstValue) / (size-1)) * 100;
			return resolution;
		}
	}else{
		return 1;
	}
}

void NetCdfConfigureDialog::createMesh()
{
	size_t* edgeT2Max = new size_t[currentVar->num_dims()];
	double originLon = 0, originLat = 0;
	double lastLon = 0, lastLat = 0;
	size_t sizeLon = 0, sizeLat = 0;
	getDimEdges(comboBoxDim1->currentIndex(), sizeLat, originLat, lastLat);
	getDimEdges(comboBoxDim2->currentIndex(), sizeLon, originLon, lastLon);

	for(int i=0; i < currentVar->num_dims(); i++) edgeT2Max[i]=1;

	// set array edges: lat x lon
	edgeT2Max[comboBoxDim1->currentIndex()]=sizeLat;
	edgeT2Max[comboBoxDim2->currentIndex()]=sizeLon;

	// set up array
	double* dimArrayT2mMax = new double[sizeLat*sizeLon];
	for(int i=0; i < (sizeLat*sizeLon); i++) dimArrayT2mMax[i]=0;

	//Time-Dimension:
	if (currentVar->num_dims() > 2)
	{
		long* newOrigin = new long[currentVar->num_dims()];
		for (int i=0; i < currentVar->num_dims(); i++) newOrigin[i]=0;
		newOrigin[comboBoxDim3->currentIndex()] = getTimeStep(); //set origin to selected time
		NcBool myStartBool = currentVar->set_cur(newOrigin);
		//Dimension 4:
		if (currentVar->num_dims() > 3) newOrigin[comboBoxDim4->currentIndex()] = getDim4(); //if there are is a 4th dimension 
		delete newOrigin;
	}
	
	NcBool arrayOfValues = currentVar->get(dimArrayT2mMax,edgeT2Max); //create Array of Values

	for (int i=0; i < (sizeLat*sizeLon); i++)
	{
		dimArrayT2mMax[i] = dimArrayT2mMax[i] - 273; // convert from kalvin to celsius
		if (dimArrayT2mMax[i] < -10000 ) dimArrayT2mMax[i] = -9999; // all values < -10000, set to "no-value"
	}
		
	std::pair <double,double> originNetCdf (originLon,originLat); // lon,lat

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

	double resolution = (doubleSpinBoxResolution->value()) * 0.01;

	currentMesh = VtkMeshConverter::convertImgToMesh(dimArrayT2mMax,originNetCdf,sizeLon,sizeLat,resolution,meshElemType,useIntensity);
	
	delete edgeT2Max;
	delete dimArrayT2mMax;
}


GridAdapter* NetCdfConfigureDialog::getMesh()
{
	return currentMesh;
}


QString NetCdfConfigureDialog::setName()
{
	std::string name;
	name.append(currentPath);
	name.erase(0,name.find_last_of("/")+1);
	name.erase(name.find_last_of("."));
	QString qstr = QString::fromStdString(name);
	return qstr;
}

std::string NetCdfConfigureDialog::getName()
{
	std::string name = (lineEditName->text()).toStdString();
	QString date = dateTimeEditDim3->date().toString(Qt::DateFormat::LocalDate);
	name.append(" - ").append(date.toStdString());
	return name;	
}




