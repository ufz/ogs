/**
 * \file mainwindow.h
 * 4/11/2009 LB Initial implementation
 *
 */

#include "mainwindow.h"

// models
#include "GEOModels.h"
#include "GeoTreeModel.h"
#include "StationTreeModel.h"
#include "MshModel.h"
#include "ConditionModel.h"

//dialogs
#include "DBConnectionDialog.h"
#include "DiagramPrefsDialog.h"
#include "GMSHPrefsDialog.h"
#include "LineEditDialog.h"
#include "ListPropertiesDialog.h"
#include "MshQualitySelectionDialog.h"
#include "SHPImportDialog.h"
#include "VtkAddFilterDialog.h"
#include "VisPrefsDialog.h"

#include "OGSRaster.h"
#include "OGSError.h"
#include "Configure.h"
#include "VtkVisPipeline.h"
#include "VtkVisPipelineItem.h"
#include "RecentFiles.h"
#include "TreeModelIterator.h"
#include "VtkGeoImageSource.h"
#include "VtkBGImageSource.h"
#include "DatabaseConnection.h"

//test
#include "rf_bc_new.h"
#include "rf_st_new.h"
#include "rf_ic_new.h"
#include "wait.h"

// FileIO includes
#include "OGSIOVer4.h"
#include "StationIO.h"
#include "PetrelInterface.h"
#include "GocadInterface.h"
#include "XMLInterface.h"
#include "GMSHInterface.h"
#include "GMSInterface.h"
#include "NetCDFInterface.h"    //YW  07.2010

//test
#include "MathIO/CRSIO.h"

// MSH
#include "msh_mesh.h"

// Qt includes
#include <QFileDialog>
#include <QMessageBox>
#include <QSettings>
#include <QDesktopWidget>

// VTK includes
#include <vtkVRMLExporter.h>
#include <vtkOBJExporter.h>

#ifdef OGS_USE_OPENSG
#include <OpenSG/OSGSceneFileHandler.h>
#include <OpenSG/OSGCoredNodePtr.h>
#include <OpenSG/OSGGroup.h>
#include "vtkOsgActor.h"
#include "OsgWidget.h"
#endif

#ifdef OGS_USE_VRPN
#include "TrackingSettingsWidget.h"
#include "VtkTrackedCamera.h"
#endif // OGS_USE_VRPN

/// FEM. 11.03.2010. WW
#include "problem.h"
Problem *aproblem = NULL;

using namespace FileIO;

MainWindow::MainWindow(QWidget *parent /* = 0*/)
: QMainWindow(parent), _db (NULL)
{
    setupUi(this);

	// Setup various models
	_geoModels = new GEOModels();
	geoTabWidget->treeView->setModel(_geoModels->getGeoModel());

	stationTabWidget->treeView->setModel(_geoModels->getStationModel());

	_meshModels = new MshModel(_project);
	mshTabWidget->treeView->setModel(_meshModels);

	_conditionModel = new ConditionModel(_project);

	// vtk visualization pipeline
#ifdef OGS_USE_OPENSG
	OsgWidget* osgWidget = new OsgWidget(this, 0, Qt::Window);
	//osgWidget->show();
	osgWidget->sceneManager()->setRoot(makeCoredNode<OSG::Group>());
	osgWidget->sceneManager()->showAll();
	_vtkVisPipeline = new VtkVisPipeline(visualizationWidget->renderer(), osgWidget->sceneManager());
#else // OGS_USE_OPENSG
	_vtkVisPipeline = new VtkVisPipeline(visualizationWidget->renderer());
#endif // OGS_USE_OPENSG

	// station model connects
	connect(stationTabWidget->treeView,
			SIGNAL(stationListExportRequested(std::string, std::string)),
			this, SLOT(exportBoreholesToGMS(std::string, std::string))); // export Stationlist to GMS
	connect(stationTabWidget->treeView,
			SIGNAL(stationListRemoved(std::string)), _geoModels,
			SLOT(removeStationVec(std::string))); // update model when stations are removed
	connect(stationTabWidget->treeView,
			SIGNAL(stationListSaved(QString, QString)), this,
			SLOT(writeStationListToFile(QString, QString))); // save Stationlist to File
	connect(_geoModels,
			SIGNAL(stationVectorRemoved(StationTreeModel*, std::string)),
			this, SLOT(updateDataViews())); // update data view when stations are removed
	connect(stationTabWidget->treeView, SIGNAL(diagramRequested(QModelIndex&)),
			this, SLOT(showDiagramPrefsDialog(QModelIndex&))); // connect treeview to diagramview

	// geo model connects
	connect(geoTabWidget->treeView, SIGNAL(listRemoved(std::string, GEOLIB::GEOTYPE)),
		_geoModels,	SLOT(removeGeometry(std::string, GEOLIB::GEOTYPE)));
	connect(geoTabWidget->treeView,	SIGNAL(saveToFileRequested(QString, QString)),
		this, SLOT(writeGeometryToFile(QString, QString))); // save geometry to file
	connect(geoTabWidget->treeView,	SIGNAL(requestLineEditDialog(const std::string&)),
		this, SLOT(showLineEditDialog(const std::string&))); // open line edit dialog
	connect(geoTabWidget->treeView,	SIGNAL(loadFEMCondFileRequested(std::string)),
		this, SLOT(loadFEMConditionsFromFile(std::string))); // add FEM Conditions
	connect(_geoModels, SIGNAL(geoDataAdded(GeoTreeModel*, std::string, GEOLIB::GEOTYPE)),
		this, SLOT(updateDataViews()));
	connect(_geoModels, SIGNAL(geoDataRemoved(GeoTreeModel*, std::string, GEOLIB::GEOTYPE)),
		this, SLOT(updateDataViews()));
	//connect(_geoModels, SIGNAL(geoDataRemoved(GeoTreeModel*, std::string, GEOLIB::GEOTYPE)),
	//	_conditionModel, SLOT(removeFEMConditions(std::string, GEOLIB::GEOTYPE)));


	// Setup connections for mesh models to GUI
	connect(mshTabWidget, SIGNAL(requestMeshRemoval(const QModelIndex&)),
			_meshModels, SLOT(removeMesh(const QModelIndex&)));
	connect(mshTabWidget, SIGNAL(qualityCheckRequested(VtkMeshSource*)),
			this, SLOT(showMshQualitySelectionDialog(VtkMeshSource*)));


	// Setup connections for condition model to GUI
	conditionTabWidget->treeView->setModel(_conditionModel);
	connect(conditionTabWidget, SIGNAL(requestConditionRemoval(const QModelIndex&)),
			_conditionModel, SLOT(removeCondition(const QModelIndex&)));

	// VisPipeline Connects
	connect(_geoModels, SIGNAL(geoDataAdded(GeoTreeModel*, std::string, GEOLIB::GEOTYPE)),
			_vtkVisPipeline, SLOT(addPipelineItem(GeoTreeModel*, std::string, GEOLIB::GEOTYPE)));
	connect(_geoModels, SIGNAL(geoDataRemoved(GeoTreeModel*, std::string, GEOLIB::GEOTYPE)),
			_vtkVisPipeline, SLOT(removeSourceItem(GeoTreeModel*, std::string, GEOLIB::GEOTYPE)));

	connect(_geoModels, SIGNAL(stationVectorAdded(StationTreeModel*, std::string)),
			_vtkVisPipeline, SLOT(addPipelineItem(StationTreeModel*, std::string)));
	connect(_geoModels, SIGNAL(stationVectorRemoved(StationTreeModel*, std::string)),
			_vtkVisPipeline, SLOT(removeSourceItem(StationTreeModel*, std::string)));

	connect(_meshModels, SIGNAL(meshAdded(MshModel*, QModelIndex)),
			_vtkVisPipeline, SLOT(addPipelineItem(MshModel*,QModelIndex)));
	connect(_meshModels, SIGNAL(meshRemoved(MshModel*, QModelIndex)),
			_vtkVisPipeline, SLOT(removeSourceItem(MshModel*, QModelIndex)));

	connect(_vtkVisPipeline, SIGNAL(vtkVisPipelineChanged()),
			visualizationWidget->vtkWidget, SLOT(update()));
	connect(_vtkVisPipeline, SIGNAL(vtkVisPipelineChanged()),
			vtkVisTabWidget->vtkVisPipelineView, SLOT(expandAll()));

	vtkVisTabWidget->vtkVisPipelineView->setModel(_vtkVisPipeline);
	connect(vtkVisTabWidget->vtkVisPipelineView,
			SIGNAL(requestRemovePipelineItem(QModelIndex)), _vtkVisPipeline,
			SLOT(removePipelineItem(QModelIndex)));
	connect(vtkVisTabWidget->vtkVisPipelineView,
			SIGNAL(requestAddPipelineFilterItem(QModelIndex)), this,
			SLOT(showAddPipelineFilterItemDialog(QModelIndex)));
	connect(vtkVisTabWidget, SIGNAL(requestViewUpdate()), visualizationWidget,
			SLOT(updateView()));

	connect(vtkVisTabWidget->vtkVisPipelineView,
			SIGNAL(actorSelected(vtkProp3D*)),
			(QObject*) (visualizationWidget->interactorStyle()),
			SLOT(highlightActor(vtkProp3D*)));

	// Propagates selected vtk object in the pipeline to the pick interactor
	connect(vtkVisTabWidget->vtkVisPipelineView,
			SIGNAL(dataObjectSelected(vtkDataObject*)),
			(QObject*) (visualizationWidget->interactorStyle()),
			SLOT(pickableDataObject(vtkDataObject*)));
	connect((QObject*) (visualizationWidget->vtkPickCallback()),
			SIGNAL(actorPicked(vtkProp3D*)),
			vtkVisTabWidget->vtkVisPipelineView, SLOT(selectItem(vtkProp3D*)));

	connect(vtkVisTabWidget->vtkVisPipelineView,
			SIGNAL(meshAdded(Mesh_Group::CFEMesh*, std::string&)),
			_meshModels, SLOT(addMesh(Mesh_Group::CFEMesh*, std::string&)));

	// Stack the data dock widgets together
	tabifyDockWidget(geoDock, mshDock);
	tabifyDockWidget(mshDock, conditionDock);
	tabifyDockWidget(conditionDock, stationDock);

	// Restore window geometry
	readSettings();

	// Get info on screens geometry(ies)
	_vtkWidget = visualizationWidget->vtkWidget;
	QDesktopWidget* desktopWidget = QApplication::desktop();
#if OGS_QT_VERSION < 46
	const unsigned int screenCount = desktopWidget->numScreens();
#else
	const unsigned int screenCount = desktopWidget->screenCount();
#endif // OGS_QT_VERSION < 46
	for (size_t i = 0; i < screenCount; ++i)
		_screenGeometries.push_back(desktopWidget->availableGeometry(i));

	// Setup import files menu
	menu_File->insertMenu(action_Exit, createImportFilesMenu());

	// Setup recent files menu
	RecentFiles* recentFiles = new RecentFiles(this, SLOT(openRecentFile()),
			"recentFileList", "OpenGeoSys-5");
	connect(this, SIGNAL(fileUsed(QString)), recentFiles,
			SLOT(setCurrentFile(QString)));
	menu_File->insertMenu(action_Exit, recentFiles->menu());

	// Setup Windows menu
	QAction* showGeoDockAction = geoDock->toggleViewAction();
	showGeoDockAction->setStatusTip(tr("Shows / hides the geometry view"));
	connect(showGeoDockAction, SIGNAL(triggered(bool)), this,
			SLOT(showGeoDockWidget(bool)));
	menuWindows->addAction(showGeoDockAction);

	QAction* showStationDockAction = stationDock->toggleViewAction();
	showStationDockAction->setStatusTip(tr("Shows / hides the station view"));
	connect(showStationDockAction, SIGNAL(triggered(bool)), this,
			SLOT(showStationDockWidget(bool)));
	menuWindows->addAction(showStationDockAction);

	QAction* showMshDockAction = mshDock->toggleViewAction();
	showMshDockAction->setStatusTip(tr("Shows / hides the mesh view"));
	connect(showMshDockAction, SIGNAL(triggered(bool)), this,
			SLOT(showMshDockWidget(bool)));
	menuWindows->addAction(showMshDockAction);

	QAction* showCondDockAction = conditionDock->toggleViewAction();
	showCondDockAction->setStatusTip(tr("Shows / hides the FEM conditions view"));
	connect(showCondDockAction, SIGNAL(triggered(bool)), this,
			SLOT(showMshDockWidget(bool)));
	menuWindows->addAction(showMshDockAction);

	QAction* showVisDockAction = vtkVisDock->toggleViewAction();
	showVisDockAction->setStatusTip(tr("Shows / hides the VTK Pipeline view"));
	connect(showVisDockAction, SIGNAL(triggered(bool)), this,
			SLOT(showVisDockWidget(bool)));
	menuWindows->addAction(showVisDockAction);

	// Presentation mode
	QMenu* presentationMenu = new QMenu();
	presentationMenu->setTitle("Presentation on");
	connect(presentationMenu, SIGNAL(aboutToShow()), this,
			SLOT(createPresentationMenu()));
	menuWindows->insertMenu(showVisDockAction, presentationMenu);

	_fileFinder.addDirectory(".");
	_fileFinder.addDirectory(std::string(SOURCEPATH).append("/FileIO"));

	#ifdef OGS_USE_VRPN
		VtkTrackedCamera* cam = static_cast<VtkTrackedCamera*>
			(visualizationWidget->renderer()->GetActiveCamera());
		_trackingSettingsWidget = new TrackingSettingsWidget(cam, visualizationWidget, Qt::Window);
	#endif // OGS_USE_VRPN

	// connects for station model
	connect(stationTabWidget->treeView,
			SIGNAL(propertiesDialogRequested(std::string)), this,
			SLOT(showPropertiesDialog(std::string)));
	//	std::cout << "size of Point: " << sizeof (GEOLIB::Point) << std::endl;
	//	std::cout << "size of CGLPoint: " << sizeof (CGLPoint) << std::endl;
	//
	//	std::cout << "size of Polyline: " << sizeof (GEOLIB::Polyline) << std::endl;
	//	std::cout << "size of CGLPolyline: " << sizeof (CGLPolyline) << std::endl;
	//
	//	std::cout << "size of GEOLIB::Surface: " << sizeof (GEOLIB::Surface) << std::endl;
	//	std::cout << "size of Surface: " << sizeof (Surface) << std::endl;
	//
	//	std::cout << "size of CCore: " << sizeof (Mesh_Group::CCore) << std::endl;
	//	std::cout << "size of CNode: " << sizeof (Mesh_Group::CNode) << std::endl;
	//	std::cout << "size of CElement: " << sizeof (Mesh_Group::CNode) << std::endl;
	//	std::cout << "size of CEdge: " << sizeof (Mesh_Group::CEdge) << std::endl;
	//	std::cout << "size of CFEMesh: " << sizeof (Mesh_Group::CFEMesh) << std::endl;
	//	std::cout << "size of Matrix: " << sizeof (Math_Group::Matrix) << std::endl;
	//
	//	std::cout << "size of vec<size_t>: " << sizeof (Math_Group::vec<size_t>) << std::endl;
	//	std::cout << "size of std::vector: " << sizeof (std::vector<size_t>) << std::endl;

	//	std::cout << "size of CSourceTerm: " << sizeof (CSourceTerm) << std::endl;
	//	std::cout << "size of CBoundaryCondition: " << sizeof (CBoundaryCondition) << std::endl;

	//	std::cout << "size of CElement: " << sizeof (FiniteElement::CElement) << std::endl;
	//	std::cout << "size of CRFProcess: " << sizeof (CRFProcess) << std::endl;
	//	std::cout << "size of CFEMesh: " << sizeof (Mesh_Group::CFEMesh) << std::endl;
}

MainWindow::~MainWindow()
{
	delete _db;
	delete _vtkVisPipeline;
	delete _meshModels;
	delete _geoModels;

#ifdef OGS_USE_VRPN
	delete _trackingSettingsWidget;
#endif // OGS_USE_VRPN
}

void MainWindow::closeEvent(QCloseEvent* event)
{
	writeSettings();
	QWidget::closeEvent(event);
}

void MainWindow::showGeoDockWidget(bool show)
{
	if (show)
		geoDock->show();
	else
		geoDock->hide();
}

void MainWindow::showStationDockWidget(bool show)
{
	if (show)
		stationDock->show();
	else
		stationDock->hide();
}

void MainWindow::showMshDockWidget(bool show)
{
	if (show)
		mshDock->show();
	else
		mshDock->hide();
}

void MainWindow::showConditionDockWidget(bool show)
{
	if (show)
		conditionDock->show();
	else
		conditionDock->hide();
}

void MainWindow::showVisDockWidget(bool show)
{
	if (show)
		vtkVisDock->show();
	else
		vtkVisDock->hide();
}

void MainWindow::open()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName = QFileDialog::getOpenFileName( this, "Select data file to open",
							settings.value("lastOpenedFileDirectory").toString(),
							"Geosys files (*.gsp *.gli *.gml *.msh *.stn);;Project files (*.gsp);;GLI files (*.gli);;MSH files (*.msh);;STN files (*.stn);;All files (* *.*)");
	if (!fileName.isEmpty()) {
		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
		loadFile(fileName);
	}
}

void MainWindow::openDatabase()
{
	if (_db == NULL) {
		_db = new DatabaseConnection(_geoModels);
		_db->dbConnect();
	}

	if (_db != NULL && _db->isConnected()) {
		_db->getListSelection();
		updateDataViews();
	}
}

void MainWindow::openDatabaseConnection()
{
	if (_db == NULL) _db = new DatabaseConnection(_geoModels);
	DBConnectionDialog* dbConn = new DBConnectionDialog();
	connect(
			dbConn,
			SIGNAL(connectionRequested(QString, QString, QString, QString, QString)),
			_db,
			SLOT(setConnection(QString, QString, QString, QString, QString)));
	dbConn->show();
}

void MainWindow::openRecentFile()
{
	QAction* action = qobject_cast<QAction*> (sender());
	if (action) loadFile(action->data().toString());
}

void MainWindow::save()
{
	QString dir_str = this->getLastUsedDir();

	QString fileName = QFileDialog::getSaveFileName(this, "Save data as", dir_str,
						"GeoSys project (*.gsp);;GeoSys4 geometry files (*.gli);;GMSH geometry files (*.geo)");

	if (!fileName.isEmpty()) {
		QFileInfo fi(fileName);

		if (fi.suffix().toLower() == "gsp") {
			std::string schemaName(_fileFinder.getPath("OpenGeoSysProject.xsd"));
			XMLInterface xml(_geoModels, schemaName);
			xml.writeProjectFile(fileName);
		/*
		} else if (fi.suffix().toLower() == "gml") {
			std::string schemaName(_fileFinder.getPath("OpenGeoSysGLI.xsd"));
			XMLInterface xml(_geoModels, schemaName);
			xml.writeGLIFile(fileName, gliName);
		*/
		} else if (fi.suffix().toLower() == "geo") {
			GMSHInterface gmsh_io(fileName.toStdString());
			std::vector<std::string> selected_geometries;
			const size_t param1(2);
			const double param2(0.3);
			const double param3(0.05);
			gmsh_io.writeAllDataToGMSHInputFile(*_geoModels,
					selected_geometries, param1, param2, param3);
		} else if (fi.suffix().toLower() == "gli") {
			//			writeGLIFileV4 (fileName.toStdString(), gliName.toStdString(), *_geoModels);
			writeAllDataToGLIFileV4(fileName.toStdString(), *_geoModels);
		}
	}
}

void MainWindow::loadFile(const QString &fileName)
{
	QFile file(fileName);
	if (!file.open(QFile::ReadOnly)) {
		QMessageBox::warning(this, tr("Application"), tr(
				"Cannot read file %1:\n%2.") .arg(fileName) .arg(
				file.errorString()));
		return;
	}

	QApplication::setOverrideCursor(Qt::WaitCursor);
	QFileInfo fi(fileName);
	std::string
			base =
					fi.absoluteDir().absoluteFilePath(fi.completeBaseName()).toStdString();
	if (fi.suffix().toLower() == "gli") {
#ifndef NDEBUG
		QTime myTimer0;
		myTimer0.start();
#endif
		//    	FileIO::readGLIFileV4 (fileName.toStdString(), _geoModels);
		readGLIFileV4(fileName.toStdString(), _geoModels);
#ifndef NDEBUG
		std::cout << myTimer0.elapsed() << " ms" << std::endl;
#endif
		//
		//#ifndef NDEBUG
		//    	QTime myTimer;
		//    	myTimer.start();
		//    	std::cout << "GEOLIB_Read_GeoLib ... " << std::flush;
		//#endif
		//    	GEOLIB_Read_GeoLib(base); //fileName.toStdString());
		//        cout << "Nr. Points: " << gli_points_vector.size() << endl;
		//		cout << "Nr. Lines: " << polyline_vector.size() << endl;
		//		cout << "Nr. Surfaces: " << surface_vector.size() << endl;
		//#ifndef NDEBUG
		//    	 std::cout << myTimer.elapsed() << " ms" << std::endl;
		//#endif
		// 		GEOCalcPointMinMaxCoordinates();
	} else if (fi.suffix().toLower() == "gsp") {
		std::string schemaName(_fileFinder.getPath("OpenGeoSysProject.xsd"));
		XMLInterface xml(_geoModels, schemaName);
		xml.readProjectFile(fileName);
	} else if (fi.suffix().toLower() == "gml") {
#ifndef NDEBUG
		QTime myTimer0;
		myTimer0.start();
#endif
		std::string schemaName(_fileFinder.getPath("OpenGeoSysGLI.xsd"));
		XMLInterface xml(_geoModels, schemaName);
		xml.readGLIFile(fileName);
#ifndef NDEBUG
		std::cout << myTimer0.elapsed() << " ms" << std::endl;
#endif
	}
	// OpenGeoSys observation station files (incl. boreholes)
	else if (fi.suffix().toLower() == "stn") {
		std::string schemaName(_fileFinder.getPath("OpenGeoSysSTN.xsd"));
		XMLInterface xml(_geoModels, schemaName);
		xml.readSTNFile(fileName);
	}
	// OpenGeoSys mesh files
	else if (fi.suffix().toLower() == "msh") {
		std::string name = fileName.toStdString();
		Mesh_Group::CFEMesh* msh = MshModel::loadMeshFromFile(name);
		if (msh)
			_meshModels->addMesh(msh, name);
		else
			OGSError::box("Failed to load a mesh file.");
	}

	// GMS borehole files
	else if (fi.suffix().toLower() == "txt") {
		std::vector<GEOLIB::Point*> *boreholes =
				new std::vector<GEOLIB::Point*>();
		std::string name = fi.baseName().toStdString();

		if (GMSInterface::readBoreholesFromGMS(boreholes,
				fileName.toStdString())) _geoModels->addStationVec(boreholes,
				name, GEOLIB::getRandomColor());
	}
	// GMS mesh files
	else if (fi.suffix().toLower() == "3dm") {
		std::string name = fileName.toStdString();
		Mesh_Group::CFEMesh* mesh = GMSInterface::readGMS3DMMesh(name);
		if (mesh) _meshModels->addMesh(mesh, name);
	}
	// goCAD files
	else if (fi.suffix().toLower() == "ts") {
#ifndef NDEBUG
		QTime myTimer;
		myTimer.start();
		std::cout << "GoCad Read ... " << std::flush;
#endif
		FileIO::GocadInterface(fileName.toStdString(), _geoModels);
#ifndef NDEBUG
		std::cout << myTimer.elapsed() << " ms" << std::endl;
#endif
	}

	// NetCDF files
	// YW  07.2010
	else if (fi.suffix().toLower() == "nc") {
#ifndef NDEBUG
		QTime myTimer;
		myTimer.start();
		std::cout << "NetCDF Read ...\n" << std::flush;
#endif
		std::string name = fileName.toStdString();
		std::vector<GEOLIB::Point*> *pnt_vec =
				new std::vector<GEOLIB::Point*>();
		/* Data dimensions. */
		size_t len_rlat, len_rlon;
		FileIO::NetCDFInterface::readNetCDFData(name, pnt_vec, _geoModels,
				len_rlat, len_rlon);
		Mesh_Group::CFEMesh* mesh = FileIO::NetCDFInterface::createMeshFromPoints(pnt_vec,
				len_rlat, len_rlon);
		//GridAdapter* grid = new GridAdapter(mesh);
		_meshModels->addMesh(mesh, name);
#ifndef NDEBUG
		std::cout << myTimer.elapsed() << " ms" << std::endl;
#endif
	}
	updateDataViews();

	emit fileUsed(fileName);
}

void MainWindow::loadPetrelFiles(const QStringList &sfc_file_names,
		const QStringList &well_path_file_names)
{
	QStringList::const_iterator it = sfc_file_names.begin();
	std::list<std::string> sfc_files;
	while (it != sfc_file_names.end()) {
		sfc_files.push_back((*it).toStdString());
		++it;
	}

	it = well_path_file_names.begin();
	std::list<std::string> well_path_files;
	while (it != well_path_file_names.end()) {
		well_path_files.push_back((*it).toStdString());
		++it;
	}

	std::string unique_str(*(sfc_files.begin()));

	PetrelInterface(sfc_files, well_path_files, unique_str, _geoModels);
}

void MainWindow::updateDataViews()
{
	visualizationWidget->showAll();
	geoTabWidget-> treeView->updateView();
	stationTabWidget-> treeView->updateView();
	mshTabWidget-> treeView->updateView();

	QApplication::restoreOverrideCursor();
}

void MainWindow::readSettings()
{
	QSettings settings("UFZ", "OpenGeoSys-5");

	restoreGeometry(settings.value("windowGeometry").toByteArray());
	restoreState(settings.value("windowState").toByteArray());
}

void MainWindow::writeSettings()
{
	QSettings settings("UFZ", "OpenGeoSys-5");

	settings.setValue("windowGeometry", saveGeometry());
	settings.setValue("windowState", saveState());
}

void MainWindow::about()
{
	QString ogsVersion = QString(OGS_VERSION);
	QMessageBox::about(this, tr("About OpenGeoSys-5"), tr(
			"Built on %1\nOGS Version: %2"). arg(
			QDate::currentDate().toString()).arg(ogsVersion));
}

QMenu* MainWindow::createImportFilesMenu()
{
	QMenu* importFiles = new QMenu("&Import Files");
	QAction* gmsFiles = importFiles->addAction("G&MS Files...");
	connect(gmsFiles, SIGNAL(triggered()), this, SLOT(importGMS()));
	QAction* gocadFiles = importFiles->addAction("&Gocad Files...");
	QAction* netcdfFiles = importFiles->addAction("&NetCDF Files...");
	connect(netcdfFiles, SIGNAL(triggered()), this, SLOT(importNetcdf()));
	connect(gocadFiles, SIGNAL(triggered()), this, SLOT(importGoCad()));
	QAction* petrelFiles = importFiles->addAction("&Petrel Files...");
	connect(petrelFiles, SIGNAL(triggered()), this, SLOT(importPetrel()));
	QAction* rasterFiles = importFiles->addAction("&Raster Files...");
	connect(rasterFiles, SIGNAL(triggered()), this, SLOT(importRaster()));
#ifdef OGS_USE_OPENSG
	QAction* rasterPolyFiles = importFiles->addAction("R&aster Files as PolyData...");
	connect(rasterPolyFiles, SIGNAL(triggered()), this, SLOT(importRasterAsPoly()));
#endif
	QAction* shapeFiles = importFiles->addAction("&Shape Files...");
	connect(shapeFiles, SIGNAL(triggered()), this, SLOT(importShape()));
	QAction* vtkFiles = importFiles->addAction("&VTK Files...");
	connect( vtkFiles, SIGNAL(triggered()), this, SLOT(importVtk()) );

	return importFiles;
}

void MainWindow::importGMS()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName = QFileDialog::getOpenFileName(this,
			"Select GMS file to import", settings.value(
					"lastOpenedFileDirectory").toString(),
			"GMS files (*.txt *.3dm)");
	if (!fileName.isEmpty()) {
		loadFile(fileName);
		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
	}
}

void MainWindow::importGoCad()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName = QFileDialog::getOpenFileName(this,
			"Select data file to import", settings.value(
					"lastOpenedFileDirectory").toString(),
			"Gocad files (*.ts);;Gocad lines (*.tline)");
	if (!fileName.isEmpty()) {
		loadFile(fileName);
		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
	}
}

void MainWindow::importRaster()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName = QFileDialog::getOpenFileName(this,
			"Select raster file to import", settings.value(
					"lastOpenedFileDirectory").toString(),
			"Raster files (*.asc *.bmp *.jpg *.png *.tif);;");
	QFileInfo fi(fileName);
	QString fileType = fi.suffix().toLower();

	if ((fileType == "asc") || (fileType == "tif") || (fileType == "png")
			|| (fileType == "jpg") || (fileType == "bmp")) {
		VtkGeoImageSource* geoImage = VtkGeoImageSource::New();
		geoImage->setImageFilename(fileName);
		_vtkVisPipeline->addPipelineItem(geoImage);

		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
	}
	else if (fileName.length() > 0) OGSError::box(
			"File extension not supported.");
}

void MainWindow::importRasterAsPoly()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName = QFileDialog::getOpenFileName(this,
			"Select raster file to import", settings.value(
					"lastOpenedFileDirectory").toString(),
			"Raster files (*.asc *.bmp *.jpg *.png *.tif);;");
	QFileInfo fi(fileName);

	if ((fi.suffix().toLower() == "asc") || (fi.suffix().toLower() == "tif")
			|| (fi.suffix().toLower() == "png") || (fi.suffix().toLower()
			== "jpg") || (fi.suffix().toLower() == "bmp")) {
		QImage raster;
		QPointF origin;
		double scalingFactor;
		OGSRaster::loadImage(fileName, raster, origin, scalingFactor, true);

		VtkBGImageSource* bg = VtkBGImageSource::New();
		bg->SetOrigin(origin.x(), origin.y());
		bg->SetCellSize(scalingFactor);
		bg->SetRaster(raster);
		bg->SetName(fileName);
		_vtkVisPipeline->addPipelineItem(bg);

		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedFileDirectory", dir.absolutePath());

	} else
		OGSError::box("File extension not supported.");
}

void MainWindow::importShape()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName = QFileDialog::getOpenFileName(this,
			"Select shape file to import", settings.value(
					"lastOpenedFileDirectory").toString(),
			"ESRI Shape files (*.shp );;");
	QFileInfo fi(fileName);

	if (fi.suffix().toLower() == "shp" || fi.suffix().toLower() == "dbf") {
		SHPImportDialog dlg((fileName.toUtf8()).constData(), _geoModels);
		dlg.exec();

		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
	}
}

void MainWindow::importPetrel()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QStringList sfc_file_names = QFileDialog::getOpenFileNames(this,
			"Select surface data file(s) to import", "", "Petrel files (*)");
	QStringList well_path_file_names = QFileDialog::getOpenFileNames(this,
			"Select well path data file(s) to import", "", "Petrel files (*)");
	if (sfc_file_names.size() != 0 || well_path_file_names.size() != 0) {
		loadPetrelFiles(sfc_file_names, well_path_file_names);
		QDir dir = QDir(sfc_file_names.at(0));
		settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
	}
}

//YW  07.2010
void MainWindow::importNetcdf()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName = QFileDialog::getOpenFileName(this,
			"Select NetCDF file to import", settings.value(
					"lastOpenedFileDirectory").toString(),
			"NetCDF files (*.nc);;");
	if (!fileName.isEmpty()) {
		loadFile(fileName);
		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
	}
}

void MainWindow::importVtk()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QStringList fileNames = QFileDialog::getOpenFileNames(this,
			"Select VTK file(s) to import", settings.value(
					"lastOpenedFileDirectory").toString(),
			"VTK files (*.vtk *.vti *.vtr *.vts *.vtp *.vtu);;");
	foreach(QString fileName, fileNames) {
		if (!fileName.isEmpty()) {
			_vtkVisPipeline->loadFromFile(fileName);
			QDir dir = QDir(fileName);
			settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
		}
	}
}

void MainWindow::showPropertiesDialog(std::string name)
{
	ListPropertiesDialog dlg(name, _geoModels);
	connect(
			&dlg,
			SIGNAL(propertyBoundariesChanged(std::string, std::vector<PropertyBounds>)),
			_geoModels,
			SLOT(filterStationVec(std::string, std::vector<PropertyBounds>)));
	dlg.exec();
}

void MainWindow::showAddPipelineFilterItemDialog(QModelIndex parentIndex)
{
	VtkAddFilterDialog dlg(_vtkVisPipeline, parentIndex);
	dlg.exec();
}

void MainWindow::loadFEMConditionsFromFile(std::string geoName)
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName = QFileDialog::getOpenFileName( this, "Select data file to open",
							settings.value("lastOpenedFileDirectory").toString(),
							"Geosys FEM condition files (*.cnd *.bc *.ic *.st);;All files (* *.*)");
	if (!fileName.isEmpty())
	{
		QFileInfo fi(fileName);
		std::vector<FEMCondition*> conditions;

		if (fi.suffix().toLower() == "cnd") {
			std::string schemaName(_fileFinder.getPath("OpenGeoSysCond.xsd"));
			XMLInterface xml(_geoModels, schemaName);
			xml.readFEMCondFile(conditions, fileName, QString::fromStdString(geoName));
			if (!conditions.empty())
				this->_conditionModel->addConditions(conditions);
		}
		else if (fi.suffix().toLower() == "bc")
		{
			QString name = fi.path() + "/";
			BCRead((name.append(fi.baseName())).toStdString(), *_geoModels, geoName);
			for (std::list<CBoundaryCondition*>::iterator it = bc_list.begin(); it != bc_list.end(); ++it)
			{
				BoundaryCondition* bc = new BoundaryCondition(*(*it), geoName);
				conditions.push_back(bc);
			}
		}
		else if (fi.suffix().toLower() == "ic")
		{
			QString name = fi.path() + "/";
			ICRead((name.append(fi.baseName())).toStdString(), *_geoModels, geoName);
			for (std::vector<CInitialCondition*>::iterator it = ic_vector.begin(); it != ic_vector.end(); ++it)
			{
				InitialCondition* ic = new InitialCondition(*(*it), geoName);
				conditions.push_back(ic);
			}
		}
		else if (fi.suffix().toLower() == "st")
		{
			QString name = fi.path() + "/";
			STRead((name.append(fi.baseName())).toStdString(), *_geoModels, geoName);
			for (std::vector<CSourceTerm*>::iterator it = st_vector.begin(); it != st_vector.end(); ++it)
			{
				SourceTerm* st = new SourceTerm(*(*it), geoName);
				conditions.push_back(st);
			}
		}

		if (!conditions.empty())
			this->_conditionModel->addConditions(conditions);
	}
}

void MainWindow::writeGeometryToFile(QString gliName, QString fileName)
{
	std::string schemaName(_fileFinder.getPath("OpenGeoSysGLI.xsd"));
	XMLInterface xml(_geoModels, schemaName);
	xml.writeGLIFile(fileName, gliName);
}

void MainWindow::writeStationListToFile(QString listName, QString fileName)
{
	std::string schemaName(_fileFinder.getPath("OpenGeoSysSTN.xsd"));
	XMLInterface xml(_geoModels, schemaName);
	xml.writeSTNFile(fileName, listName);
}

void MainWindow::exportBoreholesToGMS(std::string listName,
		std::string fileName)
{
	const std::vector<GEOLIB::Point*> *stations(_geoModels->getStationVec(
			listName));
	GMSInterface::writeBoreholesToGMS(stations, fileName);
}

void MainWindow::callGMSH(std::vector<std::string> const & selectedGeometries,
		size_t param1, double param2, double param3, double param4,
		bool delete_geo_file)
{
	if (!selectedGeometries.empty()) {
		std::cout << "Start meshing..." << std::endl;

		QString fileName("");
		QString dir_str = this->getLastUsedDir();

		if (!delete_geo_file)
			fileName = QFileDialog::getSaveFileName(this, "Save GMSH-file as",
					dir_str, "GMSH geometry files (*.geo)");
		else
			fileName = "tmp_gmsh.geo";

		if (!fileName.isEmpty())
		{
			GMSHInterface gmsh_io(fileName.toStdString());

			if (param4 == -1) { // adaptive meshing selected
				gmsh_io.writeAllDataToGMSHInputFile(*_geoModels,
						selectedGeometries, param1, param2, param3);
			} else { // homogeneous meshing selected
				gmsh_io.writeAllDataToGMSHInputFile(*_geoModels,
						selectedGeometries, param4);
			}

			if (system(NULL) != 0) { // command processor available
				std::string gmsh_command("gmsh -2 ");
				std::string fname (fileName.toStdString());
				gmsh_command += fname;
				size_t pos (fname.rfind ("."));
				if (pos != std::string::npos)
					fname = fname.substr (0, pos);
				gmsh_command += " -o " + fname + ".msh";
				system(gmsh_command.c_str());
				this->loadFile(fileName.left(fileName.length()-3).append("msh"));
			} else {
				OGSError::box("Error executing command", "Error");
			}

			if (delete_geo_file) { // delete file
				std::string remove_command ("rm ");
	#ifdef _WIN32
				remove_command = "del ";
	#endif
				remove_command += fileName.toStdString();
				std::cout << "remove command: " << remove_command << std::endl;
				system(remove_command.c_str());
			}
		}
	}
	else
	{
		OGSError::box("No geometry information selected.", "Error");
		std::cout << "No geometry information selected..." << std::endl;
	}
}

void MainWindow::showDiagramPrefsDialog(QModelIndex &index)
{
	QString listName;
	GEOLIB::Station* stn = _geoModels->getStationModel()->stationFromIndex(
			index, listName);

	if (stn->type() == GEOLIB::Station::STATION) {
		DiagramPrefsDialog prefs(stn, listName, _db);
		prefs.show();
	}
	if (stn->type() == GEOLIB::Station::BOREHOLE) OGSError::box(
			"No time series data available for borehole.");
}

void MainWindow::showLineEditDialog(const std::string &geoName)
{
	LineEditDialog lineEdit(*(_geoModels->getPolylineVecObj(geoName)));
	connect(&lineEdit, SIGNAL(connectPolylines(const std::string&, std::vector<size_t>, double, std::string, bool, bool)),
		_geoModels, SLOT(connectPolylineSegments(const std::string&, std::vector<size_t>, double, std::string, bool, bool)));
	lineEdit.exec();
}

void MainWindow::showGMSHPrefsDialog()
{
	GMSHPrefsDialog dlg(_geoModels);
	connect(
			&dlg, SIGNAL(requestMeshing(std::vector<std::string> const &, size_t, double, double, double, bool)),
			this, SLOT(callGMSH(std::vector<std::string> const &, size_t, double, double, double, bool)));
	dlg.exec();
}

void MainWindow::showMshQualitySelectionDialog(VtkMeshSource* mshSource)
{
	MshQualitySelectionDialog dlg(mshSource);
	connect(&dlg, SIGNAL(measureSelected(VtkMeshSource*, MshQualityType::type)),
			_vtkVisPipeline, SLOT(checkMeshQuality(VtkMeshSource*, MshQualityType::type)));
	dlg.exec();
}

void MainWindow::showVisalizationPrefsDialog()
{
	VisPrefsDialog dlg(_vtkVisPipeline);
	dlg.exec();
}

void MainWindow::FEMTestStart()
{
	std::string project_name ("Circle");
	std::vector<GEOLIB::Point*> *pnts (new std::vector<GEOLIB::Point*>);
	GEOLIB::Point middle_pnt (0.0, 0.0, 0.0);
	size_t resolution (36);
	GEOLIB::Polygon *polygon(createPolygonFromCircle (middle_pnt, 2.0, *pnts, resolution));
	std::vector<GEOLIB::Polyline*> *plys (new std::vector<GEOLIB::Polyline*>);
	std::cout << "polygon: " << std::endl;
	std::cout << *polygon << std::endl;
	plys->push_back (polygon);
	_geoModels->addPointVec (pnts, project_name);
	_geoModels->addPolylineVec (plys, project_name);

#ifndef NDEBUG
	std::cout << "FEM Test here ..." << std::endl;
	QSettings settings("UFZ", "OpenGeoSys-5");

	QString fileName = QFileDialog::getOpenFileName(this,
			"Select matrix file in binary compressed row storage format", settings.value(
					"lastOpenedFileDirectory").toString(),
			"binary matrix file (*.bin);;");

	std::string fname (fileName.toStdString());
	// open input stream
	std::ifstream in (fname.c_str(), std::ios::binary);

	if (in) {
		unsigned n(0), *iA(NULL), *jA(NULL);
		double *A(NULL);

		std::cout << "reading matrix ... " << std::flush;
		// read matrix
		FileIO::readCompressedStorageFmt (in, n, iA, jA, A);
		in.close ();
		std::cout << "done" << std::endl;

		// *** ToDo
		// read right hand side
		// solve system of linear equations


//		std::cout << "n : " << n << std::endl;
//		std::cout << "iA[n]: " << iA[n] << std::endl;
//		std::cout << "iA: " << std::endl;
//		for (size_t i(0); i<=n; i++) {
//			std::cout << " " << iA[i];
//		}
//		std::cout << std::endl << "jA: " << std::endl;
//		for (size_t i(0); i<iA[n]; i++) {
//			std::cout << " " << jA[i];
//		}
//		std::cout << std::endl << "A: " << std::endl;
//		for (size_t i(0); i<iA[n]; i++) {
//			std::cout << " " << A[i];
//		}
//		std::cout << std::endl;

		delete [] iA;
		delete [] jA;
		delete [] A;
	}

#else
	std::cout << "This is test functionality only..." << std::endl;
#endif
}

void MainWindow::showTrackingSettingsDialog()
{
#ifdef OGS_USE_VRPN
	_trackingSettingsWidget->show();
#else // OGS_USE_VRPN
	QMessageBox::warning(this, "Functionality not implemented",
			"Sorry but this progam was not compiled with VRPN support.");
#endif // OGS_USE_VRPN
}

void MainWindow::ShowWindow()
{
	this->show();
}

void MainWindow::HideWindow()
{
	this->hide();
}

void MainWindow::on_actionExportVTK_triggered(bool checked /*= false*/)
{
	Q_UNUSED(checked)
	QSettings settings("UFZ", "OpenGeoSys-5");
	int count = 0;
	QString filename = QFileDialog::getSaveFileName(this,
			"Export object to vtk-files", settings.value(
					"lastExportedFileDirectory").toString(),
			"VTK files (*.vtp *.vtu)");
	if (!filename.isEmpty()) {
		QDir dir = QDir(filename);
		settings.setValue("lastExportedFileDirectory", dir.absolutePath());

		std::string basename = QFileInfo(filename).path().toStdString();
		basename.append("/" + QFileInfo(filename).baseName().toStdString());
		TreeModelIterator it(_vtkVisPipeline);
		++it;
		while (*it) {
			count++;
			static_cast<VtkVisPipelineItem*> (*it)->writeToFile(basename
					+ number2str(count));
			++it;
		}
	}
}

void MainWindow::on_actionExportVRML2_triggered(bool checked /*= false*/)
{
	Q_UNUSED(checked)
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName = QFileDialog::getSaveFileName(this,
			"Save scene to VRML file", settings.value(
					"lastExportedFileDirectory").toString(),
			"VRML files (*.wrl);;");
	if (!fileName.isEmpty()) {
		QDir dir = QDir(fileName);
		settings.setValue("lastExportedFileDirectory", dir.absolutePath());

		vtkVRMLExporter* exporter = vtkVRMLExporter::New();
		exporter->SetFileName(fileName.toStdString().c_str());
		exporter->SetRenderWindow(
				visualizationWidget->vtkWidget->GetRenderWindow());
		exporter->Write();
		exporter->Delete();
	}
}

void MainWindow::on_actionExportObj_triggered(bool checked /*= false*/)
{
	Q_UNUSED(checked)
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName = QFileDialog::getSaveFileName(this,
			"Save scene to Wavefront OBJ files", settings.value(
					"lastExportedFileDirectory").toString(), ";;");
	if (!fileName.isEmpty()) {
		QDir dir = QDir(fileName);
		settings.setValue("lastExportedFileDirectory", dir.absolutePath());

		vtkOBJExporter* exporter = vtkOBJExporter::New();
		exporter->SetFilePrefix(fileName.toStdString().c_str());
		exporter->SetRenderWindow(
				visualizationWidget->vtkWidget->GetRenderWindow());
		exporter->Write();
		exporter->Delete();
	}
}

void MainWindow::on_actionExportOpenSG_triggered(bool checked /*= false*/)
{
	Q_UNUSED(checked)
#ifdef OGS_USE_OPENSG
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString filename = QFileDialog::getSaveFileName(
			this, "Export scene to OpenSG binary file", settings.value(
					"lastExportedFileDirectory").toString(), "OpenSG files (*.osb);;");
	if (!filename.isEmpty())
	{
		QDir dir = QDir(filename);
		settings.setValue("lastExportedFileDirectory", dir.absolutePath());

		TreeModelIterator it(_vtkVisPipeline);
		++it;
		OSG::NodePtr root = OSG::makeCoredNode<OSG::Group>();
		while(*it)
		{
			VtkVisPipelineItem* item = static_cast<VtkVisPipelineItem*>(*it);
			vtkOsgActor* actor = static_cast<vtkOsgActor*>(item->actor());
			actor->SetVerbose(true);
			actor->UpdateOsg();
			beginEditCP(root);
			root->addChild(actor->GetOsgRoot());
			endEditCP(root);
			actor->ClearOsg();
			++it;
		}

		OSG::SceneFileHandler::the().write(root, filename.toStdString().c_str());
	}
#else
	QMessageBox::warning(this, "Functionality not implemented",
			"Sorry but this progam was not compiled with OpenSG support.");
#endif
}

void MainWindow::createPresentationMenu()
{
	QMenu* menu = static_cast<QMenu*> (QObject::sender());
	menu->clear();
	if (!_vtkWidget->parent()) {
		QAction* action = new QAction("Quit presentation mode", menu);
		connect(action, SIGNAL(triggered()), this, SLOT(quitPresentationMode()));
		action->setShortcutContext(Qt::WidgetShortcut);
		action->setShortcut(QKeySequence(Qt::Key_Escape));
		menu->addAction(action);
	} else {
		int count = 0;
		const int currentScreen = QApplication::desktop()->screenNumber(
				visualizationWidget);
		foreach (QRect screenGeo, _screenGeometries)
			{
				Q_UNUSED(screenGeo);
				QAction* action = new QAction(
						QString("On screen %1").arg(count), menu);
				connect(action, SIGNAL(triggered()), this,
						SLOT(startPresentationMode()));
				if (count == currentScreen) action->setEnabled(false);
				menu->addAction(action);
				++count;
			}
	}
}

void MainWindow::startPresentationMode()
{
	// Save the QMainWindow state to restore when quitting presentation mode
	_windowState = this->saveState();

	// Get the screen number from the QAction which sent the signal
	QString actionText = static_cast<QAction*> (QObject::sender())->text();
	int screen = actionText.split(" ").back().toInt();

	// Move the widget to the screen and maximize it
	// Real fullscreen hides the menu
	_vtkWidget->setParent(NULL, Qt::Window);
	_vtkWidget->move(QPoint(_screenGeometries[screen].x(),
			_screenGeometries[screen].y()));
	//_vtkWidget->showFullScreen();
	_vtkWidget->showMaximized();

	// Create an action which quits the presentation mode when pressing
	// ESCAPE when the the window has focus
	QAction* action = new QAction("Quit presentation mode", this);
	connect(action, SIGNAL(triggered()), this, SLOT(quitPresentationMode()));
	action->setShortcutContext(Qt::WidgetShortcut);
	action->setShortcut(QKeySequence(Qt::Key_Escape));
	_vtkWidget->addAction(action);

	// Hide the central widget to maximize the dock widgets
	QMainWindow::centralWidget()->hide();
}

void MainWindow::quitPresentationMode()
{
	// Remove the quit action
	QAction* action = _vtkWidget->actions().back();
	_vtkWidget->removeAction(action);
	delete action;

	// Add the widget back to visualization widget
	visualizationWidget->layout()->addWidget(_vtkWidget);

	QMainWindow::centralWidget()->show();

	// Restore the previously saved QMainWindow state
	this->restoreState(_windowState);
}

QString MainWindow::getLastUsedDir()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName("");
	QStringList files = settings.value("recentFileList").toStringList();
	if (files.size() != 0)
		return QFileInfo(files[0]).absolutePath();
	else
		return QDir::homePath();
}
