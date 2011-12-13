/**
 * \file mainwindow.h
 * 4/11/2009 LB Initial implementation
 *
 */
#include "Configure.h"
#include "mainwindow.h"

// models
#include "ProcessModel.h"
#include "ElementTreeModel.h"
#include "GEOModels.h"
#include "GeoTreeModel.h"
#include "MshModel.h"
#include "StationTreeModel.h"

//dialogs
#include "DBConnectionDialog.h"
#include "DiagramPrefsDialog.h"
#include "FEMConditionSetupDialog.h"
#include "GMSHPrefsDialog.h"
#include "LineEditDialog.h"
#include "ListPropertiesDialog.h"
#include "MshQualitySelectionDialog.h"
#include "NewProcessDialog.h"
#include "SetNameDialog.h"
#include "VisPrefsDialog.h"
#include "VtkAddFilterDialog.h"

#ifdef Shapelib_FOUND
#include "SHPImportDialog.h"
#endif

#include "DatabaseConnection.h"
#include "OGSError.h"
#include "OGSRaster.h"
#include "RecentFiles.h"
#include "TreeModelIterator.h"
#include "VtkBGImageSource.h"
#include "VtkGeoImageSource.h"
#include "VtkVisPipeline.h"
#include "VtkVisPipelineItem.h"

//test
#include "BoundaryCondition.h"
#include "InitialCondition.h"
#include "MathIO/CRSIO.h"
#include "Raster.h"
#include "SourceTerm.h"
#include "rf_bc_new.h"
#include "rf_ic_new.h"
#include "rf_st_new.h"

// FileIO includes
#include "FEFLOWInterface.h"
#include "GMSInterface.h"
#include "GeoIO/Gmsh2GeoIO.h"
#include "GocadInterface.h"
#include "MeshIO/GMSHInterface.h"
#include "MeshIO/TetGenInterface.h"
#include "NetCDFInterface.h"    //YW  07.2010
#include "OGSIOVer4.h"
#include "PetrelInterface.h"
#include "StationIO.h"
#include "XmlIO/XmlCndInterface.h"
#include "XmlIO/XmlGmlInterface.h"
#include "XmlIO/XmlGspInterface.h"
#include "XmlIO/XmlStnInterface.h"

#include "StringTools.h"

// MSH
#include "msh_mesh.h"

// MSHGEOTOOLS
#include "ExtractMeshNodes.h"

// Qt includes
#include <QDesktopWidget>
#include <QFileDialog>
#include <QMessageBox>
#include <QObject>
#include <QSettings>

// VTK includes
#include <vtkOBJExporter.h>
#include <vtkRenderer.h>
#include <vtkVRMLExporter.h>

#ifdef OGS_USE_OPENSG
#include "vtkOsgConverter.h"
#include <OpenSG/OSGCoredNodePtr.h>
#include <OpenSG/OSGGroup.h>
#include <OpenSG/OSGSceneFileHandler.h>
#endif

#ifdef OGS_USE_VRPN
#include "TrackingSettingsWidget.h"
#include "VtkTrackedCamera.h"
#endif // OGS_USE_VRPN

#ifdef OGS_BUILD_INFO
#include "BuildInfo.h"
#endif // OGS_BUILD_INFO

//// test only
//#include "rf_mmp_new.h"
//#include "rf_msp_new.h"
#include "rf_mfp_new.h"

/// FEM. 11.03.2010. WW
#include "problem.h"
Problem* aproblem = NULL;

using namespace FileIO;

MainWindow::MainWindow(QWidget* parent /* = 0*/)
	: QMainWindow(parent), _db (NULL), _project()
{
	setupUi(this);

	// Setup various models
	_geoModels = new GEOModels();
	_project.setGEOObjects(_geoModels);
	_meshModels = new MshModel(_project);
	_elementModel = new ElementTreeModel();
	_processModel = new ProcessModel(_project);

	geoTabWidget->treeView->setModel(_geoModels->getGeoModel());
	stationTabWidget->treeView->setModel(_geoModels->getStationModel());
	mshTabWidget->treeView->setModel(_meshModels);
	mshTabWidget->elementView->setModel(_elementModel);
	modellingTabWidget->treeView->setModel(_processModel);

	// vtk visualization pipeline
	_vtkVisPipeline = new VtkVisPipeline(visualizationWidget->renderer());

	// station model connects
	connect(stationTabWidget->treeView, SIGNAL(stationListExportRequested(std::string, std::string)),
	        this, SLOT(exportBoreholesToGMS(std::string, std::string))); // export Stationlist to GMS
	connect(stationTabWidget->treeView, SIGNAL(stationListRemoved(std::string)), _geoModels,
	        SLOT(removeStationVec(std::string))); // update model when stations are removed
	connect(stationTabWidget->treeView, SIGNAL(stationListSaved(QString, QString)), this,
	        SLOT(writeStationListToFile(QString, QString))); // save Stationlist to File
	connect(_geoModels, SIGNAL(stationVectorRemoved(StationTreeModel *, std::string)),
	        this, SLOT(updateDataViews())); // update data view when stations are removed
	connect(stationTabWidget->treeView, SIGNAL(diagramRequested(QModelIndex &)),
	        this, SLOT(showDiagramPrefsDialog(QModelIndex &))); // connect treeview to diagramview

	// geo model connects
	connect(geoTabWidget->treeView, SIGNAL(listRemoved(std::string, GEOLIB::GEOTYPE)),
	        _geoModels, SLOT(removeGeometry(std::string, GEOLIB::GEOTYPE)));
	connect(geoTabWidget->treeView, SIGNAL(saveToFileRequested(QString, QString)),
	        this, SLOT(writeGeometryToFile(QString, QString))); // save geometry to file
	connect(geoTabWidget->treeView, SIGNAL(requestLineEditDialog(const std::string &)),
	        this, SLOT(showLineEditDialog(const std::string &))); // open line edit dialog
	connect(geoTabWidget->treeView, SIGNAL(requestNameChangeDialog(const std::string&, const GEOLIB::GEOTYPE, size_t)),
			this, SLOT(showGeoNameDialog(const std::string&, const GEOLIB::GEOTYPE, size_t)));
	connect(geoTabWidget->treeView, SIGNAL(requestCondSetupDialog(const std::string&, const GEOLIB::GEOTYPE, size_t)),
			this, SLOT(showCondSetupDialog(const std::string&, const GEOLIB::GEOTYPE, size_t)));
	connect(geoTabWidget->treeView, SIGNAL(loadFEMCondFileRequested(std::string)),
	        this, SLOT(loadFEMConditions(std::string))); // add FEM Conditions
	connect(geoTabWidget->treeView, SIGNAL(saveFEMConditionsRequested(QString, QString)),
	        this, SLOT(writeFEMConditionsToFile(QString, QString)));
	connect(_geoModels, SIGNAL(geoDataAdded(GeoTreeModel *, std::string, GEOLIB::GEOTYPE)),
	        this, SLOT(updateDataViews()));
	connect(_geoModels, SIGNAL(geoDataRemoved(GeoTreeModel *, std::string, GEOLIB::GEOTYPE)),
	        this, SLOT(updateDataViews()));
	connect(geoTabWidget->treeView, SIGNAL(geoItemSelected(const vtkPolyDataAlgorithm*, int)),
		    _vtkVisPipeline, SLOT(highlightGeoObject(const vtkPolyDataAlgorithm*, int)));
	connect(geoTabWidget->treeView, SIGNAL(removeGeoItemSelection()),
		    _vtkVisPipeline, SLOT(removeHighlightedGeoObject()));


	// Setup connections for mesh models to GUI
	connect(mshTabWidget->treeView, SIGNAL(requestMeshRemoval(const QModelIndex &)),
	        _meshModels, SLOT(removeMesh(const QModelIndex &)));
	connect(mshTabWidget->treeView, SIGNAL(requestMeshRemoval(const QModelIndex &)),
	        _elementModel, SLOT(clearView()));
	connect(mshTabWidget->treeView, SIGNAL(qualityCheckRequested(VtkMeshSource*)),
	        this, SLOT(showMshQualitySelectionDialog(VtkMeshSource*)));
	connect(mshTabWidget->treeView, SIGNAL(requestDIRECTSourceTerms(const std::vector<GEOLIB::Point*>*)),
	        this, SLOT(loadDIRECTSourceTerms(const std::vector<GEOLIB::Point*>*)));

	// Setup connections for process model to GUI
	connect(modellingTabWidget->treeView, SIGNAL(conditionsRemoved(const FiniteElement::ProcessType, const std::string&, const FEMCondition::CondType)),
	        _processModel, SLOT(removeFEMConditions(const FiniteElement::ProcessType, const std::string&, const FEMCondition::CondType)));
	connect(modellingTabWidget->treeView, SIGNAL(processRemoved(const FiniteElement::ProcessType)),
	        _processModel, SLOT(removeProcess(const FiniteElement::ProcessType)));
	connect(modellingTabWidget, SIGNAL(requestNewProcess()),
		    this, SLOT(showNewProcessDialog()));

	// VisPipeline Connects
	connect(_geoModels, SIGNAL(geoDataAdded(GeoTreeModel *, std::string, GEOLIB::GEOTYPE)),
	        _vtkVisPipeline, SLOT(addPipelineItem(GeoTreeModel *, std::string, GEOLIB::GEOTYPE)));
	connect(_geoModels, SIGNAL(geoDataRemoved(GeoTreeModel *, std::string, GEOLIB::GEOTYPE)),
	        _vtkVisPipeline, SLOT(removeSourceItem(GeoTreeModel *, std::string, GEOLIB::GEOTYPE)));

	connect(_processModel, SIGNAL(conditionAdded(ProcessModel *,  const FiniteElement::ProcessType, const FEMCondition::CondType)),
	        _vtkVisPipeline, SLOT(addPipelineItem(ProcessModel *,  const FiniteElement::ProcessType, const FEMCondition::CondType)));
	connect(_processModel, SIGNAL(conditionsRemoved(ProcessModel *, const FiniteElement::ProcessType, const FEMCondition::CondType)),
	        _vtkVisPipeline, SLOT(removeSourceItem(ProcessModel *, const FiniteElement::ProcessType, const FEMCondition::CondType)));

	connect(_geoModels, SIGNAL(stationVectorAdded(StationTreeModel *, std::string)),
	        _vtkVisPipeline, SLOT(addPipelineItem(StationTreeModel *, std::string)));
	connect(_geoModels, SIGNAL(stationVectorRemoved(StationTreeModel *, std::string)),
	        _vtkVisPipeline, SLOT(removeSourceItem(StationTreeModel *, std::string)));

	connect(_meshModels, SIGNAL(meshAdded(MshModel *, QModelIndex)),
	        _vtkVisPipeline, SLOT(addPipelineItem(MshModel *,QModelIndex)));
	connect(_meshModels, SIGNAL(meshRemoved(MshModel *, QModelIndex)),
	        _vtkVisPipeline, SLOT(removeSourceItem(MshModel *, QModelIndex)));

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
	connect((QObject*) (visualizationWidget->interactorStyle()),
	        SIGNAL(requestViewUpdate()),
	        visualizationWidget, SLOT(updateView()));

	// Propagates selected vtk object in the pipeline to the pick interactor
	connect(vtkVisTabWidget->vtkVisPipelineView,
	        SIGNAL(dataObjectSelected(vtkDataObject*)),
	        (QObject*) (visualizationWidget->interactorStyle()),
	        SLOT(pickableDataObject(vtkDataObject*)));
	connect((QObject*) (visualizationWidget->vtkPickCallback()),
	        SIGNAL(actorPicked(vtkProp3D*)),
	        vtkVisTabWidget->vtkVisPipelineView, SLOT(selectItem(vtkProp3D*)));
	connect((QObject*) (visualizationWidget->interactorStyle()),
	        SIGNAL(elementPicked(const GridAdapter *, const size_t)),
	        this->_elementModel, SLOT(setElement(const GridAdapter *, const size_t)));
	connect((QObject*) (visualizationWidget->interactorStyle()),
	        SIGNAL(elementPicked(const GridAdapter *, const size_t)),
	        mshTabWidget->elementView, SLOT(updateView()));

	connect(vtkVisTabWidget->vtkVisPipelineView,
	        SIGNAL(meshAdded(MeshLib::CFEMesh *, std::string &)),
	        _meshModels, SLOT(addMesh(MeshLib::CFEMesh *, std::string &)));

	// Stack the data dock widgets together
	tabifyDockWidget(geoDock, mshDock);
	tabifyDockWidget(mshDock, modellingDock);
	tabifyDockWidget(modellingDock, stationDock);

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
		_screenGeometries.push_back(desktopWidget->availableGeometry((int)i));

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

	QAction* showModellingDockAction = modellingDock->toggleViewAction();
	showModellingDockAction->setStatusTip(tr("Shows / hides the Process view"));
	connect(showModellingDockAction, SIGNAL(triggered(bool)), this,
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
#endif     // OGS_USE_VRPN

	// connects for station model
	connect(stationTabWidget->treeView,
	        SIGNAL(propertiesDialogRequested(std::string)), this,
	        SLOT(showPropertiesDialog(std::string)));

	_visPrefsDialog = new VisPrefsDialog(_vtkVisPipeline, visualizationWidget);

	//	std::cout << "size of Point: " << sizeof (GEOLIB::Point) << std::endl;
	//	std::cout << "size of CGLPoint: " << sizeof (CGLPoint) << std::endl;
	//
	//	std::cout << "size of Polyline: " << sizeof (GEOLIB::Polyline) << std::endl;
	//	std::cout << "size of CGLPolyline: " << sizeof (CGLPolyline) << std::endl;
	//
	//	std::cout << "size of GEOLIB::Surface: " << sizeof (GEOLIB::Surface) << std::endl;
	//	std::cout << "size of Surface: " << sizeof (Surface) << std::endl;
	//
	//	std::cout << "size of CCore: " << sizeof (MeshLib::CCore) << std::endl;
	//	std::cout << "size of CNode: " << sizeof (MeshLib::CNode) << std::endl;
	//	std::cout << "size of CElement: " << sizeof (MeshLib::CNode) << std::endl;
	//	std::cout << "size of CEdge: " << sizeof (MeshLib::CEdge) << std::endl;
	//	std::cout << "size of CFEMesh: " << sizeof (MeshLib::CFEMesh) << std::endl;
	//	std::cout << "size of Matrix: " << sizeof (Math_Group::Matrix) << std::endl;
	//
	//	std::cout << "size of vec<size_t>: " << sizeof (Math_Group::vec<size_t>) << std::endl;
	//	std::cout << "size of std::vector: " << sizeof (std::vector<size_t>) << std::endl;

	//	std::cout << "size of CSourceTerm: " << sizeof (CSourceTerm) << std::endl;
	//	std::cout << "size of CBoundaryCondition: " << sizeof (CBoundaryCondition) << std::endl;

	//	std::cout << "size of CElement: " << sizeof (FiniteElement::CElement) << std::endl;
//		std::cout << "size of CMediumProperties: " << sizeof(CMediumProperties) << std::endl;
//		std::cout << "size of CSolidProperties: " << sizeof(SolidProp::CSolidProperties) << std::endl;
		std::cout << "size of CFluidProperties: " << sizeof(CFluidProperties) << std::endl;
	//	std::cout << "size of CRFProcess: " << sizeof (CRFProcess) << std::endl;
	//	std::cout << "size of CFEMesh: " << sizeof (MeshLib::CFEMesh) << std::endl;
}

MainWindow::~MainWindow()
{
	delete _db;
	delete _vtkVisPipeline;
	delete _meshModels;
	delete _processModel;
	//delete _visPrefsDialog;
	//delete _geoModels;

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
		modellingDock->show();
	else
		modellingDock->hide();
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
	QString fileName = QFileDialog::getOpenFileName( this, "Select data file to open",settings.value("lastOpenedFileDirectory").toString(),
	                                                 "Geosys files (*.gsp *.gli *.gml *.msh *.stn);;Project files (*.gsp);;GeoSys FEM Conditions (*.cnd *.bc *.ic *.st);;GLI files (*.gli);;MSH files (*.msh);;STN files (*.stn);;All files (* *.*)");
	if (!fileName.isEmpty())
	{
		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
		loadFile(fileName);
	}
}

void MainWindow::openDatabase()
{
	if (_db == NULL)
	{
		_db = new DatabaseConnection(_geoModels);
		_db->dbConnect();
	}

	if (_db != NULL && _db->isConnected())
	{
		_db->getListSelection();
		updateDataViews();
	}
}

void MainWindow::openDatabaseConnection()
{
	if (_db == NULL)
		_db = new DatabaseConnection(_geoModels);
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
	if (action)
		loadFile(action->data().toString());
}

void MainWindow::save()
{
	QString dir_str = this->getLastUsedDir();

	QString fileName = QFileDialog::getSaveFileName(
	        this,
	        "Save data as",
	        dir_str,
	        "GeoSys project (*.gsp);;GeoSys4 geometry files (*.gli);;GMSH geometry files (*.geo)");

	if (!fileName.isEmpty())
	{
		QFileInfo fi(fileName);

		if (fi.suffix().toLower() == "gsp")
		{
			std::string schemaName(_fileFinder.getPath("OpenGeoSysProject.xsd"));
			XmlGspInterface xml(&_project, schemaName);
			xml.writeFile(fileName);
		}
		else if (fi.suffix().toLower() == "geo")
		{
			// it works like this (none of it is particularily fast or optimised or anything):
			// 1. merge all geometries that are currently loaded, all of these will be integrated into the mesh
			// 2. if "useStationsAsConstraints"-parameter is true, GMSH-Interface will also integrate all stations that are currently loaded
			//    if "useSteinerPoints"-parameter is true, additional points will be inserted in large areas without information
			// 3. after the geo-file is created the merged geometry is deleted again as it is no longer needed
			GMSHInterface gmsh_io(fileName.toStdString());
			std::vector<std::string> names;
			this->_project.getGEOObjects()->getGeometryNames(names);
			std::string merge_name("MergedGeometry");
			_geoModels->mergeGeometries (names, merge_name);
			gmsh_io.writeGMSHInputFile(merge_name,
			                           *(this->_project.getGEOObjects()), true, true);
			this->_project.getGEOObjects()->removeSurfaceVec(merge_name);
			this->_project.getGEOObjects()->removePolylineVec(merge_name);
			this->_project.getGEOObjects()->removePointVec(merge_name);
		}
		else if (fi.suffix().toLower() == "gli")
			//			writeGLIFileV4 (fileName.toStdString(), gliName.toStdString(), *_geoModels);
			writeAllDataToGLIFileV4(fileName.toStdString(), *_geoModels);
	}
}

void MainWindow::loadFile(const QString &fileName)
{
	QFile file(fileName);
	if (!file.open(QFile::ReadOnly))
	{
		QMessageBox::warning(this, tr("Application"), tr(
		                             "Cannot read file %1:\n%2.").arg(fileName).arg(
		                             file.errorString()));
		return;
	}

	QApplication::setOverrideCursor(Qt::WaitCursor);
	QFileInfo fi(fileName);
	std::string base =
	        fi.absoluteDir().absoluteFilePath(fi.completeBaseName()).toStdString();
	if (fi.suffix().toLower() == "gli")
	{
#ifndef NDEBUG
		QTime myTimer0;
		myTimer0.start();
#endif
		//      FileIO::readGLIFileV4 (fileName.toStdString(), _geoModels);
		std::string unique_name;
		readGLIFileV4(fileName.toStdString(), _geoModels, unique_name);
#ifndef NDEBUG
		std::cout << myTimer0.elapsed() << " ms" << std::endl;
#endif
		//
		//#ifndef NDEBUG
		//      QTime myTimer;
		//      myTimer.start();
		//      std::cout << "GEOLIB_Read_GeoLib ... " << std::flush;
		//#endif
		//      GEOLIB_Read_GeoLib(base); //fileName.toStdString());
		//        cout << "Nr. Points: " << gli_points_vector.size() << endl;
		//		cout << "Nr. Lines: " << polyline_vector.size() << endl;
		//		cout << "Nr. Surfaces: " << surface_vector.size() << endl;
		//#ifndef NDEBUG
		//       std::cout << myTimer.elapsed() << " ms" << std::endl;
		//#endif
		//              GEOCalcPointMinMaxCoordinates();
	}
	else if (fi.suffix().toLower() == "gsp")
	{
		std::string schemaName(_fileFinder.getPath("OpenGeoSysProject.xsd"));
		XmlGspInterface xml(&_project, schemaName);
		xml.readFile(fileName);
		std::cout << "Adding missing meshes to GUI..." << std::endl;
		_meshModels->updateModel();
	}
	else if (fi.suffix().toLower() == "gml")
	{
#ifndef NDEBUG
		QTime myTimer0;
		myTimer0.start();
#endif
		std::string schemaName(_fileFinder.getPath("OpenGeoSysGLI.xsd"));
		XmlGmlInterface xml(&_project, schemaName);
		xml.readFile(fileName);
#ifndef NDEBUG
		std::cout << myTimer0.elapsed() << " ms" << std::endl;
#endif
	}
	// OpenGeoSys observation station files (incl. boreholes)
	else if (fi.suffix().toLower() == "stn")
	{
		std::string schemaName(_fileFinder.getPath("OpenGeoSysSTN.xsd"));
		XmlStnInterface xml(&_project, schemaName);
		xml.readFile(fileName);
	}
	// OpenGeoSys mesh files
	else if (fi.suffix().toLower() == "msh")
	{
		std::string name = fileName.toStdString();
		MeshLib::CFEMesh* msh = FileIO::OGSMeshIO::loadMeshFromFile(name);
		if (msh)
			_meshModels->addMesh(msh, name);
		else
			OGSError::box("Failed to load a mesh file.");
	}
	else if ((fi.suffix().toLower() == "cnd") ||
		     (fi.suffix().toLower() == "bc") ||
			 (fi.suffix().toLower() == "ic") ||
			 (fi.suffix().toLower() == "st"))
	{
		this->loadFEMConditionsFromFile(fileName);
	}

	// GMS borehole files
	else if (fi.suffix().toLower() == "txt")
	{
		std::vector<GEOLIB::Point*>* boreholes =
		        new std::vector<GEOLIB::Point*>();
		std::string name = fi.baseName().toStdString();

		if (GMSInterface::readBoreholesFromGMS(boreholes, fileName.toStdString()))
			_geoModels->addStationVec(boreholes, name, GEOLIB::getRandomColor());
		else
			OGSError::box("Error reading GMS file.");
	}
	// GMS mesh files
	else if (fi.suffix().toLower() == "3dm")
	{
		std::string name = fileName.toStdString();
		MeshLib::CFEMesh* mesh = GMSInterface::readGMS3DMMesh(name);
		if (mesh)
			_meshModels->addMesh(mesh, name);
	}
	// goCAD files
	else if (fi.suffix().toLower() == "ts")
	{
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
	else if (fi.suffix().toLower() == "nc")
	{
#ifndef NDEBUG
		QTime myTimer;
		myTimer.start();
		std::cout << "NetCDF Read ...\n" << std::flush;
#endif
		std::string name = fileName.toStdString();
		std::vector<GEOLIB::Point*>* pnt_vec =
		        new std::vector<GEOLIB::Point*>();
		/* Data dimensions. */
		size_t len_rlat, len_rlon;
		FileIO::NetCDFInterface::readNetCDFData(name, pnt_vec, _geoModels,
		                                        len_rlat, len_rlon);
		MeshLib::CFEMesh* mesh = FileIO::NetCDFInterface::createMeshFromPoints(pnt_vec,
		                                                                       len_rlat,
		                                                                       len_rlon);
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
	while (it != sfc_file_names.end())
	{
		sfc_files.push_back((*it).toStdString());
		++it;
	}

	it = well_path_file_names.begin();
	std::list<std::string> well_path_files;
	while (it != well_path_file_names.end())
	{
		well_path_files.push_back((*it).toStdString());
		++it;
	}

	std::string unique_str(*(sfc_files.begin()));

	PetrelInterface(sfc_files, well_path_files, unique_str, _geoModels);
}

void MainWindow::updateDataViews()
{
	visualizationWidget->updateViewOnLoad();
	geoTabWidget->treeView->updateView();
	stationTabWidget->treeView->updateView();
	mshTabWidget->treeView->updateView();

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
	
	QString about = tr("Built on %1\nOGS Version: %2\n\n").arg(
		QDate::currentDate().toString(Qt::ISODate)).arg(ogsVersion);
#ifdef OGS_BUILD_INFO
#ifdef SVN_REVISION
	about.append(QString("Svn commit: %1\n").arg(SVN_REVISION));
#endif
#ifdef GIT_COMMIT_INFO
	QString gitCommit = QString(GIT_COMMIT_INFO);
	about.append(QString("Git commit: %1\n").arg(gitCommit.mid(7)));
#endif // GIT_COMMIT_INFO
#ifdef GIT_BRANCH_INFO
	QString gitBranch = QString(GIT_BRANCH_INFO);
	about.append(QString("Git branch: %1\n").arg(gitBranch.mid(2)));
#endif // GIT_BRANCH_INFO
#endif // OGS_BUILD_INFO
	QMessageBox::about(this, "About OpenGeoSys-5", about);
}

QMenu* MainWindow::createImportFilesMenu()
{
	QMenu* importFiles = new QMenu("&Import Files");
	QAction* feflowFiles = importFiles->addAction("&FEFLOW Files...");
	connect(feflowFiles, SIGNAL(triggered()), this, SLOT(importFeflow()));
	QAction* gmsFiles = importFiles->addAction("G&MS Files...");
	connect(gmsFiles, SIGNAL(triggered()), this, SLOT(importGMS()));
	QAction* gocadFiles = importFiles->addAction("&Gocad Files...");
	connect(gocadFiles, SIGNAL(triggered()), this, SLOT(importGoCad()));
	QAction* netcdfFiles = importFiles->addAction("&NetCDF Files...");
	connect(netcdfFiles, SIGNAL(triggered()), this, SLOT(importNetcdf()));
	QAction* petrelFiles = importFiles->addAction("&Petrel Files...");
	connect(petrelFiles, SIGNAL(triggered()), this, SLOT(importPetrel()));
	QAction* rasterFiles = importFiles->addAction("&Raster Files...");
	connect(rasterFiles, SIGNAL(triggered()), this, SLOT(importRaster()));
#ifdef OGS_USE_OPENSG
	QAction* rasterPolyFiles = importFiles->addAction("R&aster Files as PolyData...");
	connect(rasterPolyFiles, SIGNAL(triggered()), this, SLOT(importRasterAsPoly()));
#endif
#ifdef Shapelib_FOUND
	QAction* shapeFiles = importFiles->addAction("&Shape Files...");
	connect(shapeFiles, SIGNAL(triggered()), this, SLOT(importShape()));
#endif
	QAction* tetgenFiles = importFiles->addAction("&TetGen Files...");
	connect( tetgenFiles, SIGNAL(triggered()), this, SLOT(importTetGen()) );
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
	if (!fileName.isEmpty())
	{
		loadFile(fileName);
		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
	}
}

void MainWindow::importGoCad()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName = QFileDialog::getOpenFileName(this,
	                                                "Select data file to import",
	                                                settings.value(
	                                                        "lastOpenedFileDirectory").toString(),
	                                                "Gocad files (*.ts);;Gocad lines (*.tline)");
	if (!fileName.isEmpty())
	{
		loadFile(fileName);
		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
	}
}

void MainWindow::importRaster()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
#ifdef libgeotiff_FOUND
	QString geotiffExtension(" *.tif");
#else
	QString geotiffExtension("");
#endif
	QString fileName = QFileDialog::getOpenFileName(this, "Select raster file to import",
					settings.value("lastOpenedFileDirectory").toString(), QString(
									"Raster files (*.asc *.bmp *.jpg *.png%1);;") .arg(geotiffExtension));

	if (!fileName.isEmpty())
	{
		VtkGeoImageSource* geoImage = VtkGeoImageSource::New();
		geoImage->setImageFilename(fileName);
		_vtkVisPipeline->addPipelineItem(geoImage);

		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
	}
}

void MainWindow::importRasterAsPoly()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
#ifdef libgeotiff_FOUND
	QString geotiffExtension(" *.tif");
#else
	QString geotiffExtension("");
#endif
	QString fileName = QFileDialog::getOpenFileName(this, "Select raster file to import",
					settings.value("lastOpenedFileDirectory").toString(), QString(
									"Raster files (*.asc *.bmp *.jpg *.png%1);;") .arg(
									geotiffExtension));

	if (!fileName.isEmpty())
	{
		QImage raster;
		QPointF origin;
		double scalingFactor;
		OGSRaster::loadImage(fileName, raster, origin, scalingFactor, false);
		//OGSRaster::loadImage(fileName, raster, origin, scalingFactor, true, true);

		VtkBGImageSource* bg = VtkBGImageSource::New();
		bg->SetOrigin(origin.x(), origin.y());
		bg->SetCellSize(scalingFactor);
		bg->SetRaster(raster);
		bg->SetName(fileName);
		_vtkVisPipeline->addPipelineItem(bg);

		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
	}
}

#ifdef Shapelib_FOUND
void MainWindow::importShape()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName = QFileDialog::getOpenFileName(this, "Select shape file to import",
					settings.value("lastOpenedFileDirectory").toString(),
	                                                "ESRI Shape files (*.shp );;");
	QFileInfo fi(fileName);

	if (fi.suffix().toLower() == "shp" || fi.suffix().toLower() == "dbf")
	{
		SHPImportDialog dlg((fileName.toUtf8()).constData(), _geoModels);
		dlg.exec();

		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
	}
}
#endif

void MainWindow::importPetrel()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QStringList sfc_file_names = QFileDialog::getOpenFileNames(
	        this, "Select surface data file(s) to import", "", "Petrel files (*)");
	QStringList well_path_file_names = QFileDialog::getOpenFileNames(
	        this, "Select well path data file(s) to import", "", "Petrel files (*)");
	if (sfc_file_names.size() != 0 || well_path_file_names.size() != 0)
	{
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
	                                                "Select NetCDF file to import",
	                                                settings.value(
	                                                        "lastOpenedFileDirectory").toString(),
	                                                "NetCDF files (*.nc);;");
	if (!fileName.isEmpty())
	{
		loadFile(fileName);
		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
	}
}

void MainWindow::importTetGen()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName = QFileDialog::getOpenFileName(this,
	                                                "Select TetGen file to import",
	                                                settings.value(
	                                                        "lastOpenedFileDirectory").toString(),
	                                                "TetGen files (*.nc);;");
	if (!fileName.isEmpty())
	{
		loadFile(fileName);
		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
	}}

void MainWindow::importVtk()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QStringList fileNames = QFileDialog::getOpenFileNames(this,
	                                                      "Select VTK file(s) to import",
	                                                      settings.value(
	                                                              "lastOpenedFileDirectory").
	                                                      toString(),
	                                                      "VTK files (*.vtk *.vti *.vtr *.vts *.vtp *.vtu);;");
	foreach(QString fileName, fileNames) {
		if (!fileName.isEmpty())
		{
			_vtkVisPipeline->loadFromFile(fileName);
			QDir dir = QDir(fileName);
			settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
		}
	}
}

void MainWindow::importFeflow()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName = QFileDialog::getOpenFileName(this,
	                                                "Select FEFLOW file(s) to import",
	                                                settings.value(
	                                                        "lastOpenedFileDirectory").toString(),
	                                                "FEFLOW files (*.fem);;");
	if (!fileName.isEmpty())
	{
		FEFLOWInterface feflowIO(_geoModels);
		MeshLib::CFEMesh* msh = feflowIO.readFEFLOWModelFile(fileName.toStdString());
		if (msh)
		{
			std::string str = fileName.toStdString();
			_meshModels->addMesh(msh, str);
			QDir dir = QDir(fileName);
			settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
			//_geoModels->modified("Feflow");
			updateDataViews();
		}
		else
			OGSError::box("Failed to load a FEFLOW file.");
	}
	emit fileUsed(fileName);
}

void MainWindow::showPropertiesDialog(std::string const& name)
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


void MainWindow::loadFEMConditions(std::string geoName)
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName = QFileDialog::getOpenFileName( this, "Select data file to open",
														settings.value(
																"lastOpenedFileDirectory").
														toString(),
														"Geosys FEM condition files (*.cnd *.bc *.ic *.st);;All files (* *.*)");
	QDir dir = QDir(fileName);
	settings.setValue("lastOpenedFileDirectory", dir.absolutePath());

	if (!fileName.isEmpty())
		this->loadFEMConditionsFromFile(fileName, geoName);
}

void MainWindow::loadFEMConditionsFromFile(const QString &fileName, std::string geoName)
{
	std::vector<FEMCondition*> conditions;
	QFileInfo fi(fileName);
	if (fi.suffix().toLower() == "cnd")
	{
		std::string schemaName(_fileFinder.getPath("OpenGeoSysCond.xsd"));
		XmlCndInterface xml(&_project, schemaName);
		xml.readFile(conditions, fileName);
	}
	else
	{
		if (geoName.empty())
		{
			// assume that geoName is identical to filename of the currently loaded file (but with *.gli-extension)
			QFileInfo fi(fileName);
			geoName = fi.fileName().toStdString();
			geoName = geoName.substr(0, geoName.find_last_of(".")).append(".gli");
		}
		if (fi.suffix().toLower() == "bc")
		{
			QString name = fi.path() + "/";
			BCRead((name.append(fi.baseName())).toStdString(), *_geoModels, geoName);
			for (std::list<CBoundaryCondition*>::iterator it = bc_list.begin();
					it != bc_list.end(); ++it)
			{
				BoundaryCondition* bc = new BoundaryCondition(*(*it), geoName);
				conditions.push_back(bc);
			}
		}
		else if (fi.suffix().toLower() == "ic")
		{
			QString name = fi.path() + "/";
			ICRead((name.append(fi.baseName())).toStdString(), *_geoModels, geoName);
			for (std::vector<CInitialCondition*>::iterator it = ic_vector.begin();
					it != ic_vector.end(); ++it)
			{
				InitialCondition* ic = new InitialCondition(*(*it), geoName);
				conditions.push_back(ic);
			}
		}
		else if (fi.suffix().toLower() == "st")
		{
			QString name = fi.path() + "/";
			STRead((name.append(fi.baseName())).toStdString(), *_geoModels, geoName);
			for (std::vector<CSourceTerm*>::iterator it = st_vector.begin();
					it != st_vector.end(); ++it)
			{
				SourceTerm* st = new SourceTerm(*(*it), geoName);
				conditions.push_back(st);
			}
		}
	}
	if (!conditions.empty())
	{
		this->_processModel->addConditions(conditions);

		for (std::list<CBoundaryCondition*>::iterator it = bc_list.begin();
			    it != bc_list.end(); ++it)
			delete *it;
		bc_list.clear();
		for (size_t i = 0; i < ic_vector.size(); i++)
			delete ic_vector[i];
		ic_vector.clear();
		for (size_t i = 0; i < st_vector.size(); i++)
			delete st_vector[i];
		st_vector.clear();
	}
}

void MainWindow::writeFEMConditionsToFile(QString geoName, QString fileName)
{
	std::string schemaName(_fileFinder.getPath("OpenGeoSysCond.xsd"));
	XmlCndInterface xml(&_project, schemaName);
	xml.writeFile(fileName, geoName);
}

void MainWindow::writeGeometryToFile(QString gliName, QString fileName)
{
	std::string schemaName(_fileFinder.getPath("OpenGeoSysGLI.xsd"));
	XmlGmlInterface xml(&_project, schemaName);
	xml.writeFile(fileName, gliName);
}

void MainWindow::writeStationListToFile(QString listName, QString fileName)
{
	std::string schemaName(_fileFinder.getPath("OpenGeoSysSTN.xsd"));
	XmlStnInterface xml(&_project, schemaName);
	xml.writeFile(fileName, listName);
}

void MainWindow::exportBoreholesToGMS(std::string listName,
                                      std::string fileName)
{
	const std::vector<GEOLIB::Point*>* stations(_geoModels->getStationVec(
	                                                    listName));
	GMSInterface::writeBoreholesToGMS(stations, fileName);
}

void MainWindow::callGMSH(std::vector<std::string> const & selectedGeometries,
                          size_t param1, double param2, double param3, double param4,
                          bool delete_geo_file)
{
	if (!selectedGeometries.empty())
	{
		std::cout << "Start meshing..." << std::endl;

		QString fileName("");
		QString dir_str = this->getLastUsedDir();

		if (!delete_geo_file)
			fileName = QFileDialog::getSaveFileName(this,
			                                        "Save GMSH-file as",
			                                        dir_str,
			                                        "GMSH geometry files (*.geo)");
		else
			fileName = "tmp_gmsh.geo";

		if (!fileName.isEmpty())
		{
			GMSHInterface gmsh_io(fileName.toStdString());

			if (param4 == -1) // adaptive meshing selected
				gmsh_io.writeAllDataToGMSHInputFile(*_geoModels,
				                                    selectedGeometries,
				                                    param1,
				                                    param2,
				                                    param3);
			else // homogeneous meshing selected
				gmsh_io.writeAllDataToGMSHInputFile(*_geoModels,
				                                    selectedGeometries, param4);

			if (system(NULL) != 0) // command processor available
			{
				std::string gmsh_command("gmsh -2 ");
				std::string fname (fileName.toStdString());
				gmsh_command += fname;
				size_t pos (fname.rfind ("."));
				if (pos != std::string::npos)
					fname = fname.substr (0, pos);
				gmsh_command += " -o " + fname + ".msh";
				system(gmsh_command.c_str());
				this->loadFile(fileName.left(fileName.length() - 3).append("msh"));
			}
			else
				OGSError::box(
				        "Error executing command gmsh - no command processor available",
				        "Error");

			if (delete_geo_file) // delete file
			{
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

	if (stn->type() == GEOLIB::Station::STATION)
	{
		DiagramPrefsDialog* prefs = new DiagramPrefsDialog(stn, listName, _db);
		prefs->setAttribute(Qt::WA_DeleteOnClose);
		prefs->show();
	}
	if (stn->type() == GEOLIB::Station::BOREHOLE)
		OGSError::box(
		        "No time series data available for borehole.");
}

void MainWindow::showDiagramPrefsDialog()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName = QFileDialog::getOpenFileName( this, "Select data file to open",
	                                                 settings.value(
	                                                         "lastOpenedFileDirectory").
	                                                 toString(),
	                                                 "Text files (*.txt);;All files (* *.*)");
	if (!fileName.isEmpty())
	{
		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
		DiagramPrefsDialog* prefs = new DiagramPrefsDialog(fileName);
		prefs->setAttribute(Qt::WA_DeleteOnClose);
		prefs->show();
	}
}

void MainWindow::showGeoNameDialog(const std::string &geometry_name, const GEOLIB::GEOTYPE object_type, size_t id)
{
	std::string old_name = this->_geoModels->getElementNameByID(geometry_name, object_type, id);
	SetNameDialog dlg(geometry_name, GEOLIB::convertGeoTypeToString(object_type), id, old_name);
	connect(&dlg, SIGNAL(requestNameChange(const std::string&, const GEOLIB::GEOTYPE, size_t, std::string)),
		this->_geoModels, SLOT(addNameForElement(const std::string&, const GEOLIB::GEOTYPE, size_t, std::string)));
	dlg.exec();

	static_cast<GeoTreeModel*>(this->geoTabWidget->treeView->model())->setNameForItem(geometry_name, object_type,
		id,	this->_geoModels->getElementNameByID(geometry_name, object_type, id));
}

void MainWindow::showCondSetupDialog(const std::string &geometry_name, const GEOLIB::GEOTYPE object_type, size_t id)
{
	std::string geo_name = this->_geoModels->getElementNameByID(geometry_name, object_type, id);
	if (geo_name.empty())
	{
		this->showGeoNameDialog(geometry_name, object_type, id);
		geo_name = this->_geoModels->getElementNameByID(geometry_name, object_type, id);
	}
	// Object should now have a name ... if not, cancel the setup process
	if (geo_name.empty())
		OGSError::box("FEM Condition Setup cancelled.");
	else
	{
		FEMConditionSetupDialog dlg(geometry_name, object_type, geo_name, this->_geoModels->getGEOObject(geometry_name, object_type, geo_name));
		connect(&dlg, SIGNAL(addFEMCondition(FEMCondition*)), this->_processModel, SLOT(addCondition(FEMCondition*)));
		dlg.exec();
	}
}

void MainWindow::showNewProcessDialog()
{
	NewProcessDialog dlg;
	connect(&dlg , SIGNAL(addProcess(ProcessInfo*)),
	        _processModel, SLOT(addProcess(ProcessInfo*)));
	dlg.exec();
}

void MainWindow::showLineEditDialog(const std::string &geoName)
{
	LineEditDialog lineEdit(*(_geoModels->getPolylineVecObj(geoName)));
	connect(&lineEdit, SIGNAL(connectPolylines(const std::string &, std::vector<size_t>, double, std::string, bool, bool)),
	        _geoModels, SLOT(connectPolylineSegments(const std::string &, std::vector<size_t>, double, std::string, bool, bool)));
	lineEdit.exec();
}

void MainWindow::showGMSHPrefsDialog()
{
	GMSHPrefsDialog dlg(_geoModels);
	connect(&dlg, SIGNAL(requestMeshing(std::vector<std::string> const &, size_t, double, double, double, bool)),
	        this, SLOT(callGMSH(std::vector<std::string> const &, size_t, double, double, double, bool)));
	dlg.exec();
}

void MainWindow::showMshQualitySelectionDialog(VtkMeshSource* mshSource)
{
	MshQualitySelectionDialog dlg(mshSource);
	connect(&dlg, SIGNAL(measureSelected(VtkMeshSource *, MshQualityType::type)),
	        _vtkVisPipeline, SLOT(checkMeshQuality(VtkMeshSource *, MshQualityType::type)));
	dlg.exec();
}

void MainWindow::showVisalizationPrefsDialog()
{
	_visPrefsDialog->show();
}

void MainWindow::FEMTestStart()
{
	//FEMConditionSetupDialog* dlg = new FEMConditionSetupDialog(this->_project);
	//dlg->exec();

	// *** begin test TetGen read mesh
	const std::string path ("/home/fischeth/Desktop/data/Ketzin/PSglobal/Tom/MSH/");
	std::string mesh_name ("ClosedSurface");
	std::string fname_nodes(path + mesh_name + ".1.node");
	std::string fname_elements(path + mesh_name + ".1.ele");

	FileIO::TetGenInterface tetgen;
	MeshLib::CFEMesh* mesh (tetgen.readTetGenMesh (fname_nodes, fname_elements));

	if (mesh)
		_meshModels->addMesh(mesh, mesh_name);
	else
		OGSError::box("Failed to load TetGen mesh file.");	std::cout << "This is test functionality only..." << std::endl;
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
	                                                "Export object to vtk-files",
	                                                settings.value(
	                                                        "lastExportedFileDirectory").
	                                                toString(),
	                                                "VTK files (*.vtp *.vtu)");
	if (!filename.isEmpty())
	{
		QDir dir = QDir(filename);
		settings.setValue("lastExportedFileDirectory", dir.absolutePath());

		std::string basename = QFileInfo(filename).path().toStdString();
		basename.append("/" + QFileInfo(filename).baseName().toStdString());
		TreeModelIterator it(_vtkVisPipeline);
		++it;
		while (*it)
		{
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
	                                                        "lastExportedFileDirectory").
	                                                toString(),
	                                                "VRML files (*.wrl);;");
	if (!fileName.isEmpty())
	{
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
	                                                "Save scene to Wavefront OBJ files",
	                                                settings.value(
	                                                        "lastExportedFileDirectory").
	                                                toString(),
	                                                ";;");
	if (!fileName.isEmpty())
	{
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
			vtkOsgConverter converter(static_cast<vtkActor*>(item->actor()));
			if(converter.WriteAnActor())
			{
				beginEditCP(root);
				root->addChild(converter.GetOsgNode());
				endEditCP(root);
			}
			++it;
		}

		OSG::SceneFileHandler::the().write(root, filename.toStdString().c_str());
	}
#else // ifdef OGS_USE_OPENSG
	QMessageBox::warning(this, "Functionality not implemented",
	                     "Sorry but this progam was not compiled with OpenSG support.");
#endif
}

void MainWindow::createPresentationMenu()
{
	QMenu* menu = static_cast<QMenu*> (QObject::sender());
	menu->clear();
	if (!_vtkWidget->parent())
	{
		QAction* action = new QAction("Quit presentation mode", menu);
		connect(action, SIGNAL(triggered()), this, SLOT(quitPresentationMode()));
		action->setShortcutContext(Qt::WidgetShortcut);
		action->setShortcut(QKeySequence(Qt::Key_Escape));
		menu->addAction(action);
	}
	else
	{
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
			if (count == currentScreen)
				action->setEnabled(false);
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

void MainWindow::loadDIRECTSourceTerms(const std::vector<GEOLIB::Point*>* points)
{
	std::string geo_name("mshNodes");

	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName = QFileDialog::getOpenFileName( this, "Select data file to open",
	                                                 settings.value(
	                                                         "lastOpenedFileDirectory").
	                                                 toString(),
	                                                 "Geosys FEM condition files (*.st);;All files (* *.*)");
	QFileInfo fi(fileName);
	QString name = fi.path() + "/";

	if (!fileName.isEmpty())
	{
		// create new geometry points vector by copying mesh nodes vector
		std::vector<GEOLIB::Point*>* new_points = new std::vector<GEOLIB::Point*>;
		std::map<std::string, size_t>* name_pnt_id_map = new std::map<std::string, size_t>;

		for (size_t i = 0; i < points->size(); i++)
		{
			GEOLIB::Point* pnt = new GEOLIB::Point((*(*points)[i])[0],
			                                       (*(*points)[i])[1],
			                                       (*(*points)[i])[2]);
			new_points->push_back(pnt);
			std::stringstream out;
			out << i;
			name_pnt_id_map->insert(std::pair<std::string, size_t>(out.str(), i));
		}
		this->_geoModels->addPointVec(new_points, geo_name, name_pnt_id_map);

		STRead((name.append(fi.baseName())).toStdString(), *_geoModels, geo_name);
		std::vector<FEMCondition*> conditions = SourceTerm::createDirectSourceTerms(st_vector,geo_name);

		// add boundary conditions to model
		if (!conditions.empty())
		{
			this->_processModel->addConditions(conditions);
			for (size_t i = 0; i < st_vector.size(); i++)
				delete st_vector[i];
			st_vector.clear();
		}
	}
}

