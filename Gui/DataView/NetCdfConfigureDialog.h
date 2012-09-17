//file NetCDFConfigureDialog.cpp
//CH Initial implementation

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
 * \brief The dialog which presents the user options to load NetCDF-files.
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
	void setVariableSelect();
	void setDimensionSelect();
	void getDimEdges(int dimId,std::size_t &size, double &firstValue, double &lastValue);
	void getDaysTime(double minSince, QTime &time, int &days);
	long convertDateToMinutes(QDateTime initialDateTime,QDate selectedDate, QTime selectedTime);
	void createDataObject();
	int valueWithMaxDim();
	int getTimeStep();
	int getDim4();
	double getResolution();
	QString setName();
	void reverseNorthSouth(double* data, std::size_t width, std::size_t height);

	NcFile *_currentFile;
	NcVar *_currentVar;
	QDateTime _currentInitialDateTime;
	MeshLib::Mesh* _currentMesh;
	VtkGeoImageSource* _currentRaster;
	std::string _currentPath;
};

#endif //NETCDFCONFIGUREDIALOG_H

