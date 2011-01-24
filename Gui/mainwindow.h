/**
* \file mainwindow.h
* 4/11/2009 LB Initial implementation
*
*/

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "ui_mainwindow.h"
#include "FileFinder.h"
#include "ProjectData.h"

class GEOModels;
class MshModel;
class StationTreeModel;
class ConditionModel;
class VtkVisPipeline;
class DatabaseConnection;

#ifdef OGS_USE_VRPN
	class TrackingSettingsWidget;
#endif // OGS_USE_VRPN

class MainWindow : public QMainWindow, public Ui_MainWindowClass
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();

	void ShowWindow();
	void HideWindow();

protected:
	void closeEvent( QCloseEvent* event );

protected slots:
	void showPntDockWidget( bool show );
	void showLineDockWidget( bool show );
	void showSurfaceDockWidget( bool show );
	void showStationDockWidget( bool show );
	void showMshDockWidget( bool show );
	void showVisDockWidget( bool show );

    void open();
	void save();
	void callGMSH(std::vector<std::string> const & selectedGeometries, size_t param1, double param2, double param3, double param4);
	void exportBoreholesToGMS(std::string listName, std::string fileName);
    void importGMS();
	void importGoCad();
	void importRaster();
	void importRasterAsPoly();
	void importShape();
	void importPetrel();
	void importNetcdf();     //YW  07.2010
	void importVtk();
	void openDatabase();
	void openDatabaseConnection();
	void openRecentFile();
	void about();
	void showDiagramPrefsDialog(QModelIndex &index);
	void showGMSHPrefsDialog();
	void showPropertiesDialog(std::string name);
	void showVisalizationPrefsDialog();
	void showTrackingSettingsDialog();
	void updateDataViews();
	void writeStationListToFile(QString listName, QString fileName);
	void showAddPipelineFilterItemDialog(QModelIndex parentIndex);

	void on_actionExportVTK_triggered(bool checked = false);
	void on_actionExportVRML2_triggered(bool checked = false);
	void on_actionExportObj_triggered(bool checked = false);
	void on_actionExportOpenSG_triggered(bool checked = false);

	void createPresentationMenu();
	void startPresentationMode();
	void quitPresentationMode();

private:
	QMenu* createImportFilesMenu();
    void loadFile(const QString &fileName);
    void loadPetrelFiles(const QStringList &sfc_file_names, const QStringList &well_path_file_names);

	void readSettings();
	void writeSettings();

    QString curFile;

	DatabaseConnection* _db;
	FileFinder _fileFinder;
	GEOModels* _geoModels;
	MshModel* _meshModels;
	ConditionModel* _conditionModel;
	ProjectData _project;
	VtkVisPipeline* _vtkVisPipeline;
	QList<QRect> _screenGeometries;
	QWidget* _vtkWidget;
	QByteArray _windowState;
	#ifdef OGS_USE_VRPN
		TrackingSettingsWidget* _trackingSettingsWidget;
	#endif // OGS_USE_VRPN

signals:
	void fileUsed( QString filename );
};

class StartQt4
{
public:
	StartQt4()
	{
		int i = 0;
		QApplication* qapp = new QApplication(i, NULL);
		qapp->exec();
	}
};

#endif // MAINWINDOW_H
