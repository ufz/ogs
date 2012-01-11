//file NetCDFConfigureDialog.cpp
//CH Initial implementation

#ifndef NETCDFCONFIGUREDIALOG_H
#define NETCDFCONFIGUREDIALOG_H

#include <vtknetcdf\netcdf.h>
#include <vtknetcdf\netcdfcpp.h>
#include <QDialog>
#include "ui_NetCdfConfigure.h"
#include "GEOObjects.h"
#include "msh_mesh.h"
#include "..\Qt\VtkVis\VtkMeshConverter.h"
#include <QSettings>
#include <QMessageBox>

class NetCdfConfigureDialog : public QDialog, private Ui_NetCdfConfigure
{
	Q_OBJECT

public:
	NetCdfConfigureDialog(char* fileName, QDialog* parent = 0);
	~NetCdfConfigureDialog(void);
	GridAdapter* getMesh();
	std::string getName();
private:
	void setVariableSelect();
	void setDimensionSelect();
	void getDimEdges(int dimId,size_t &size, double &firstValue, double &lastValue);
	void getDaysTime(double minSince, QTime &time, int &days);
	long convertDateToMinutes(QDateTime initialDateTime,QDate selectedDate, QTime selectedTime);
	void createMesh();
	int valueWithMaxDim();
	int getTimeStep();
	int getDim4();
	int getResolution();
	NcFile *currentFile;
	NcVar *currentVar;
	QDateTime currentInitialDateTime;
	GridAdapter* currentMesh;
	char* currentPath;
	
private slots:
	void accept();
	void reject();
	void on_comboBoxVariable_currentIndexChanged(int id);
	void on_comboBoxDim1_currentIndexChanged(int id);
	void on_comboBoxDim2_currentIndexChanged(int id);
	void on_comboBoxDim3_currentIndexChanged(int id);
	void on_comboBoxDim4_currentIndexChanged(int id);
	
};

#endif //NETCDFCONFIGUREDIALOG_H
