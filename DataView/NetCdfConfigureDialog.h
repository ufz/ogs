//file NetCDFConfigureDialog.cpp
//CH Initial implementation

#ifndef NETCDFCONFIGUREDIALOG_H
#define NETCDFCONFIGUREDIALOG_H

#ifdef VTK_NETCDF_FOUND
#include <vtknetcdf/netcdfcpp.h>
#else
#include <netcdfcpp.h>
#endif
#include <QDialog>
#include "ui_NetCdfConfigure.h"

class GridAdapter;
class VtkGeoImageSource;

class NetCdfConfigureDialog : public QDialog, private Ui_NetCdfConfigure
{
	Q_OBJECT

public:
	NetCdfConfigureDialog(const std::string &fileName, QDialog* parent = 0);
	~NetCdfConfigureDialog(void);
	GridAdapter* getMesh() { return _currentMesh; };
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
	void setVariableSelect();
	void setDimensionSelect();
	void getDimEdges(int dimId,size_t &size, double &firstValue, double &lastValue);
	void getDaysTime(double minSince, QTime &time, int &days);
	long convertDateToMinutes(QDateTime initialDateTime,QDate selectedDate, QTime selectedTime);
	void createDataObject();
	int valueWithMaxDim();
	int getTimeStep();
	int getDim4();
	double getResolution();
	QString setName();
	void reverseNorthSouth(double* data, size_t width, size_t height);
	
	NcFile *_currentFile;
	NcVar *_currentVar;
	QDateTime _currentInitialDateTime;
	GridAdapter* _currentMesh;
	VtkGeoImageSource* _currentRaster;
	std::string _currentPath;
};

#endif //NETCDFCONFIGUREDIALOG_H

