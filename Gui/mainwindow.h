/**
 * \file mainwindow.h
 * 4/11/2009 LB Initial implementation
 *
 */

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "FileFinder.h"
#include "ProjectData.h"
#include "ui_mainwindow.h"

class GEOModels;
class MshModel;
class ElementTreeModel;
class StationTreeModel;
class ProcessModel;
class VtkVisPipeline;
class DatabaseConnection;
class VisPrefsDialog;

#ifdef OGS_USE_VRPN
class TrackingSettingsWidget;
#endif // OGS_USE_VRPN

/**
 * Main program window for the graphical user interface of OpenGeoSys.
 */
class MainWindow : public QMainWindow, public Ui_MainWindowClass
{
	Q_OBJECT

public:
	MainWindow(QWidget* parent = 0);
	~MainWindow();

	void ShowWindow();
	void HideWindow();

protected:
	void closeEvent( QCloseEvent* event );

protected slots:
	void showGeoDockWidget( bool show );
	void showStationDockWidget( bool show );
	void showMshDockWidget( bool show );
	void showConditionDockWidget( bool show );
	void showVisDockWidget( bool show );

	/// Function calls for opening files.
	void open();
	/// Function calls for saving files.
	void save();
	/// Function calls for generating GMSH files from the GUI
	void callGMSH(std::vector<std::string> & selectedGeometries,
	              size_t param1,
	              double param2,
	              double param3,
	              double param4,
	              bool delete_geo_file);
	/// Function calls for GMS export.
	void exportBoreholesToGMS(std::string listName, std::string fileName);
	/// Testing functionality for connection to FEM lib
	void FEMTestStart();
	void importGMS();
	void importGoCad();
	void importRaster();
	void importRasterAsPoly();
#ifdef Shapelib_FOUND
	void importShape();
#endif
	void importPetrel();
	void importNetcdf(); //YW  07.2010
	void importVtk();
	void importFeflow();
	void importTetGen();
	void loadFEMConditions(std::string geoName);
	void openDatabase();
	void openDatabaseConnection();
	void openRecentFile();
	void about();
	void showAddPipelineFilterItemDialog(QModelIndex parentIndex);
	void showConditionWriterDialog();
	/// Call dialog for creating or modifying FEM conditions.
	void showCondSetupDialog(const std::string &geometry_name, const GEOLIB::GEOTYPE object_type, size_t id, bool on_points = false);
	/// Allows setting the name for a geometric object
	void showGeoNameDialog(const std::string &geometry_name, const GEOLIB::GEOTYPE object_type, size_t id);
	/// Calls the diagram prefs dialog from the Tools menu.
	void showDiagramPrefsDialog();
	/// Calls the diagram prefs dialog from the station list (i.e. for a specific station).
	void showDiagramPrefsDialog(QModelIndex &index);
	void showLineEditDialog(const std::string &geoName);
	void showGMSHPrefsDialog();
	void showMshQualitySelectionDialog(VtkMeshSource* mshSource);
	void showNewProcessDialog();
	void showPropertiesDialog(std::string const& name);
	void showVisalizationPrefsDialog();
	void showTrackingSettingsDialog();
	void updateDataViews();
	void addFEMConditions(const std::vector<FEMCondition*> conditions);
	void writeFEMConditionsToFile(const QString &geoName, const FEMCondition::CondType type, const QString &fileName);
	void writeGeometryToFile(QString listName, QString fileName);
	void writeStationListToFile(QString listName, QString fileName);

	void on_actionExportVTK_triggered(bool checked = false);
	void on_actionExportVRML2_triggered(bool checked = false);
	void on_actionExportObj_triggered(bool checked = false);
	void on_actionExportOpenSG_triggered(bool checked = false);

	void createPresentationMenu();
	void startPresentationMode();
	void quitPresentationMode();

	void loadDIRECTSourceTermsFromASCII(const std::string name, const std::vector<GEOLIB::Point*>* points); //TODO put this in a better place

private:
	QMenu* createImportFilesMenu();
	void loadFile(const QString &fileName);
	void loadFEMConditionsFromFile(const QString &fileName, std::string geoName = "");
	void loadPetrelFiles(const QStringList &sfc_file_names,
	                     const QStringList &well_path_file_names);

	void readSettings();
	void writeSettings();
	QString getLastUsedDir();

	QString curFile;

	DatabaseConnection* _db;
	FileFinder _fileFinder;
	GEOModels* _geoModels;
	MshModel* _meshModels;
	ElementTreeModel* _elementModel;
	ProcessModel* _processModel;
	ProjectData _project;
	VtkVisPipeline* _vtkVisPipeline;
	QList<QRect> _screenGeometries;
	QWidget* _vtkWidget;
	QByteArray _windowState;
	QMenu* _import_files_menu;

#ifdef OGS_USE_VRPN
	TrackingSettingsWidget* _trackingSettingsWidget;
#endif     // OGS_USE_VRPN
	VisPrefsDialog* _visPrefsDialog;

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
