/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file mainwindow.h
 *
 * Created on 2009-11-04 by Lars Bilke
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
#include "CondFromRasterDialog.h"
#include "ConditionWriterDialog.h"
#include "DiagramPrefsDialog.h"
#include "FEMConditionSetupDialog.h"
//TODO6 #include "OGSFileConverter.h"
#include "GMSHPrefsDialog.h"
#include "LineEditDialog.h"
#include "ListPropertiesDialog.h"
#include "MshQualitySelectionDialog.h"
#include "NetCdfConfigureDialog.h"
#include "NewProcessDialog.h"
#include "SetNameDialog.h"
#include "VisPrefsDialog.h"
#include "VtkAddFilterDialog.h"

#ifdef Shapelib_FOUND
#include "SHPImportDialog.h"
#endif

#include "OGSError.h"
#include "VtkRaster.h"
#include "RecentFiles.h"
#include "TreeModelIterator.h"
#include "VtkBGImageSource.h"
#include "VtkGeoImageSource.h"
#include "VtkVisPipeline.h"
#include "VtkVisPipelineItem.h"

// FEM Conditions
#include "BoundaryCondition.h"
#include "InitialCondition.h"
#include "SourceTerm.h"

// FileIO includes
// TODO6 #include "FEFLOWInterface.h"
#include "GMSInterface.h"
#include "Gmsh2GeoIO.h"
#include "Legacy/MeshIO.h"
#include "Legacy/OGSIOVer4.h"
#include "MeshIO/GMSHInterface.h"
#include "MeshIO/TetGenInterface.h"
#include "PetrelInterface.h"
#include "XmlIO/VTKInterface.h"
#include "XmlIO/XmlCndInterface.h"
#include "XmlIO/XmlGmlInterface.h"
#include "XmlIO/XmlGspInterface.h"
#include "XmlIO/XmlStnInterface.h"

#include "StringTools.h"

// MeshLib
#include "Mesh.h"
#include "Node.h"
#include "MshEditor.h"

// Qt includes
#include <QDesktopWidget>
#include <QFileDialog>
#include <QMessageBox>
#include <QObject>
#include <QSettings>
#include <QSignalMapper>

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


using namespace FileIO;

MainWindow::MainWindow(QWidget* parent /* = 0*/)
	: QMainWindow(parent), _project(), _import_files_menu(NULL)
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
	connect(geoTabWidget->treeView, SIGNAL(listRemoved(std::string, GeoLib::GEOTYPE)),
	        _geoModels, SLOT(removeGeometry(std::string, GeoLib::GEOTYPE)));
	connect(geoTabWidget->treeView, SIGNAL(saveToFileRequested(QString, QString)),
	        this, SLOT(writeGeometryToFile(QString, QString))); // save geometry to file
	connect(geoTabWidget->treeView, SIGNAL(requestLineEditDialog(const std::string &)),
	        this, SLOT(showLineEditDialog(const std::string &))); // open line edit dialog
	connect(geoTabWidget->treeView, SIGNAL(requestNameChangeDialog(const std::string&, const GeoLib::GEOTYPE, std::size_t)),
			this, SLOT(showGeoNameDialog(const std::string&, const GeoLib::GEOTYPE, std::size_t)));
	connect(geoTabWidget->treeView, SIGNAL(requestCondSetupDialog(const std::string&, const GeoLib::GEOTYPE, std::size_t, bool)),
			this, SLOT(showCondSetupDialog(const std::string&, const GeoLib::GEOTYPE, std::size_t, bool)));
	connect(geoTabWidget->treeView, SIGNAL(loadFEMCondFileRequested(std::string)),
	        this, SLOT(loadFEMConditions(std::string))); // add FEM Conditions
	//connect(geoTabWidget->treeView, SIGNAL(saveFEMConditionsRequested(QString, QString)),
	//        this, SLOT(writeFEMConditionsToFile(QString, QString)));
	connect(_geoModels, SIGNAL(geoDataAdded(GeoTreeModel *, std::string, GeoLib::GEOTYPE)),
	        this, SLOT(updateDataViews()));
	connect(_geoModels, SIGNAL(geoDataRemoved(GeoTreeModel *, std::string, GeoLib::GEOTYPE)),
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
	connect(mshTabWidget->treeView, SIGNAL(requestCondSetupDialog(const std::string&, const GeoLib::GEOTYPE, std::size_t, bool)),
			this, SLOT(showCondSetupDialog(const std::string&, const GeoLib::GEOTYPE, std::size_t, bool)));

	// Setup connections for process model to GUI
	connect(modellingTabWidget->treeView, SIGNAL(conditionsRemoved(const FiniteElement::ProcessType, const std::string&, const FEMCondition::CondType)),
	        _processModel, SLOT(removeFEMConditions(const FiniteElement::ProcessType, const std::string&, const FEMCondition::CondType)));
	connect(modellingTabWidget->treeView, SIGNAL(processRemoved(const FiniteElement::ProcessType)),
	        _processModel, SLOT(removeProcess(const FiniteElement::ProcessType)));
	connect(modellingTabWidget, SIGNAL(requestNewProcess()),
		    this, SLOT(showNewProcessDialog()));
	connect(modellingTabWidget->treeView, SIGNAL(saveConditionsRequested()),
			this, SLOT(showConditionWriterDialog()));

	// VisPipeline Connects
	connect(_geoModels, SIGNAL(geoDataAdded(GeoTreeModel *, std::string, GeoLib::GEOTYPE)),
	        _vtkVisPipeline, SLOT(addPipelineItem(GeoTreeModel *, std::string, GeoLib::GEOTYPE)));
	connect(_geoModels, SIGNAL(geoDataRemoved(GeoTreeModel *, std::string, GeoLib::GEOTYPE)),
	        _vtkVisPipeline, SLOT(removeSourceItem(GeoTreeModel *, std::string, GeoLib::GEOTYPE)));

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
	        SIGNAL(elementPicked(const MeshLib::Mesh *, const std::size_t)),
	        this->_elementModel, SLOT(setElement(const MeshLib::Mesh *, const std::size_t)));
	connect((QObject*) (visualizationWidget->interactorStyle()),
	        SIGNAL(elementPicked(const MeshLib::Mesh *, const std::size_t)),
	        mshTabWidget->elementView, SLOT(updateView()));

	connect(vtkVisTabWidget->vtkVisPipelineView, SIGNAL(meshAdded(MeshLib::Mesh*)),
	        _meshModels, SLOT(addMesh(MeshLib::Mesh*)));

	// Stack the data dock widgets together
	tabifyDockWidget(geoDock, mshDock);
	tabifyDockWidget(mshDock, modellingDock);
	tabifyDockWidget(modellingDock, stationDock);

	// Restore window geometry
	readSettings();

	// Get info on screens geometry(ies)
	_vtkWidget = visualizationWidget->vtkWidget;
	QDesktopWidget* desktopWidget = QApplication::desktop();
	const unsigned int screenCount = desktopWidget->screenCount();
	for (size_t i = 0; i < screenCount; ++i)
		_screenGeometries.push_back(desktopWidget->availableGeometry((int)i));

	// Setup import files menu
	_import_files_menu = createImportFilesMenu();
	menu_File->insertMenu(action_Exit, _import_files_menu);

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
}

MainWindow::~MainWindow()
{
	delete _signal_mapper;
	delete _import_files_menu;
	delete _vtkVisPipeline;
	delete _meshModels;
	delete _processModel;

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

void MainWindow::open(int file_type)
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	ImportFileType::type t = static_cast<ImportFileType::type>(file_type);
	QString type_str = QString::fromStdString((ImportFileType::convertImportFileTypeToString(t)));
	QString fileName = QFileDialog::getOpenFileName(this,
	                                                "Select " + type_str + " file to import",
													settings.value("lastOpenedFileDirectory").toString(),
	                                                QString::fromStdString(ImportFileType::getFileSuffixString(t)));
	if (!fileName.isEmpty())
	{
		loadFile(t, fileName);
		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
	}

}

void MainWindow::openRecentFile()
{
	QAction* action = qobject_cast<QAction*> (sender());
	if (action)
		loadFile(ImportFileType::OGS, action->data().toString());
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
			xml.writeToFile(fileName.toStdString());
		}
		else if (fi.suffix().toLower() == "geo")
		{
			// it works like this (none of it is particularily fast or optimised or anything):
			// 1. merge all geometries that are currently loaded, all of these will be integrated into the mesh
			// 2. if "useStationsAsConstraints"-parameter is true, GMSH-Interface will also integrate all stations that are currently loaded
			//    if "useSteinerPoints"-parameter is true, additional points will be inserted in large areas without information
			// 3. after the geo-file is created the merged geometry is deleted again as it is no longer needed
			std::vector<std::string> names;
			this->_project.getGEOObjects()->getGeometryNames(names);
			std::string merge_name("MergedGeometry");
			_geoModels->mergeGeometries (names, merge_name);
			names.clear();
			names.push_back(merge_name);

			double param1(0.5); // mesh density scaling on normal points
			double param2(0.05); // mesh density scaling on station points
			size_t param3(2); // points per leaf
			GMSHInterface gmsh_io(*(this->_project.getGEOObjects()), true, FileIO::GMSH::AdaptiveMeshDensity, param1, param2, param3, names);
			gmsh_io.writeToFile(fileName.toStdString());

			this->_project.getGEOObjects()->removeSurfaceVec(merge_name);
			this->_project.getGEOObjects()->removePolylineVec(merge_name);
			this->_project.getGEOObjects()->removePointVec(merge_name);
		}
	}
}

void MainWindow::loadFile(ImportFileType::type t, const QString &fileName)
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
	std::string base = fi.absoluteDir().absoluteFilePath(fi.completeBaseName()).toStdString();

	if (t == ImportFileType::OGS)
	{
		if (fi.suffix().toLower() == "gli")
		{
			std::string unique_name;
			std::vector<std::string> errors;
			if (! readGLIFileV4(fileName.toStdString(), _geoModels, unique_name, errors)) {
				for (size_t k(0); k<errors.size(); k++)
					OGSError::box(QString::fromStdString(errors[k]));
			}
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
			std::string schemaName(_fileFinder.getPath("OpenGeoSysGLI.xsd"));
			XmlGmlInterface xml(&_project, schemaName);
			xml.readFile(fileName);
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
#ifndef NDEBUG
			QTime myTimer0;
			myTimer0.start();
#endif
			FileIO::MeshIO meshIO;
			std::string name = fileName.toStdString();
			MeshLib::Mesh* msh = meshIO.loadMeshFromFile(name);
			if (msh)
			{
				_meshModels->addMesh(msh);
#ifndef NDEBUG
				std::cout << "Total mesh loading time: " << myTimer0.elapsed() << " ms" << std::endl;
#endif
			}
			else
				OGSError::box("Failed to load a mesh file.");
		}
		else if (fi.suffix().toLower() == "vtu")
		{
			MeshLib::Mesh* msh = FileIO::VTKInterface::readVTUFile(fileName.toStdString());
			if (msh)
				_meshModels->addMesh(msh);
			else
				OGSError::box("Failed to load a mesh file.");
		}
		else if (fi.suffix().toLower() == "cnd")
		{
			this->loadFEMConditionsFromFile(fileName);
		}
	}
	else if (t == ImportFileType::FEFLOW)
	{
		OGSError::box("Interface not yet integrated");
		/* TODO6
		FEFLOWInterface feflowIO(_geoModels);
		MeshLib::Mesh* msh = feflowIO.readFEFLOWModelFile(fileName.toStdString());
		if (msh)
		{
			std::string name = fileName.toStdString();
			msh->setName(name);
			_meshModels->addMesh(msh);
			QDir dir = QDir(fileName);
			settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
			updateDataViews();
		}
		else
			OGSError::box("Failed to load a FEFLOW file.");
		*/
	}
	else if (t == ImportFileType::GMS)
	{
		if (fi.suffix().toLower() == "txt") // GMS borehole files
		{
			std::vector<GeoLib::Point*>* boreholes = new std::vector<GeoLib::Point*>();
			std::string name = fi.baseName().toStdString();

			if (GMSInterface::readBoreholesFromGMS(boreholes, fileName.toStdString()))
				_geoModels->addStationVec(boreholes, name);
			else
				OGSError::box("Error reading GMS file.");
		}
		else if (fi.suffix().toLower() == "3dm") // GMS mesh files
		{
			std::string name = fileName.toStdString();
			MeshLib::Mesh* mesh = GMSInterface::readGMS3DMMesh(name);
			if (mesh)
				_meshModels->addMesh(mesh);
		}
	}
	else if (t == ImportFileType::GMSH)
	{
		std::string msh_name (fileName.toStdString());
		if (FileIO::GMSHInterface::isGMSHMeshFile (msh_name))
		{
			MeshLib::Mesh* mesh = FileIO::GMSHInterface::readGMSHMesh(msh_name);
			if (mesh)
				_meshModels->addMesh(mesh);
			return;
		}
	}
	else if (t == ImportFileType::NETCDF)	// CH  01.2012
	{
		MeshLib::Mesh* mesh (NULL);

		NetCdfConfigureDialog dlg(fileName.toStdString());
		dlg.exec();
		if (dlg.getMesh() != NULL)
		{
			mesh = dlg.getMesh();
			mesh->setName(dlg.getName());
			_meshModels->addMesh(mesh);
		}
		if (dlg.getRaster() != NULL)
		{
			_vtkVisPipeline->addPipelineItem(dlg.getRaster());
		}
	}
	else if (t == ImportFileType::RASTER)
	{
		VtkGeoImageSource* geoImage = VtkGeoImageSource::New();
		geoImage->readImage(fileName);
		_vtkVisPipeline->addPipelineItem(geoImage);
	}
	else if (t == ImportFileType::POLYRASTER)
	{
		QImage raster;
		double origin[2];
		double cellSize;
		vtkImageAlgorithm* imageAlgorithm = VtkRaster::loadImage(fileName.toStdString(), origin[0], origin[1], cellSize);
		VtkBGImageSource* bg = VtkBGImageSource::New();
		bg->SetRaster(imageAlgorithm, origin[0], origin[1], cellSize);
		bg->SetName(fileName);
		_vtkVisPipeline->addPipelineItem(bg);
		//QDir dir = QDir(fileName);
		//settings.setValue("lastOpenedRasterFileDirectory", dir.absolutePath());
	}
#ifdef Shapelib_FOUND
	else if (t == ImportFileType::SHAPE)
	{
		SHPImportDialog dlg((fileName.toUtf8()).constData(), _geoModels);
		dlg.exec();
		//QDir dir = QDir(fileName);
		//settings.setValue("lastOpenedShapeFileDirectory", dir.absolutePath());
	}
#endif
	else if (t == ImportFileType::TETGEN)
	{
		QSettings settings("UFZ", "OpenGeoSys-5");
		QString element_fname = QFileDialog::getOpenFileName(this, "Select TetGen element file",
						                                     settings.value("lastOpenedTetgenFileDirectory").toString(),
						                                     "TetGen element files (*.ele);;");

		if (!fileName.isEmpty() && !element_fname.isEmpty()) {
			FileIO::TetGenInterface tetgen;
			MeshLib::Mesh* msh (tetgen.readTetGenMesh(fileName.toStdString(), element_fname.toStdString()));
			if (msh) {
				std::string name(fileName.toStdString());
				msh->setName(name);
				_meshModels->addMesh(msh);
			} else
				OGSError::box("Failed to load a TetGen mesh.");
			settings.setValue("lastOpenedTetgenFileDirectory", QDir(fileName).absolutePath());
		}
	}
	else if (t == ImportFileType::VTK)
	{
		_vtkVisPipeline->loadFromFile(fileName);
		//QDir dir = QDir(fileName);
		//settings.setValue("lastOpenedVtkFileDirectory", dir.absolutePath());
	}
	updateDataViews();

	if (t == ImportFileType::OGS) emit fileUsed(fileName);
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
	_signal_mapper = new QSignalMapper(this);
	QMenu* importFiles = new QMenu("&Import Files");
/* TODO6
	QAction* feflowFiles = importFiles->addAction("&FEFLOW Files...");
	connect(feflowFiles, SIGNAL(triggered()), this, SLOT(importFeflow()));
*/
	QAction* gmsFiles = importFiles->addAction("G&MS Files...");
	connect(gmsFiles, SIGNAL(triggered()), _signal_mapper, SLOT(map()));
	_signal_mapper->setMapping(gmsFiles, ImportFileType::GMS);
	QAction* gmshFiles = importFiles->addAction("&GMSH Files...");
	connect(gmshFiles, SIGNAL(triggered()), _signal_mapper, SLOT(map()));
	_signal_mapper->setMapping(gmshFiles, ImportFileType::GMSH);
	QAction* netcdfFiles = importFiles->addAction("&NetCDF Files...");
	connect(netcdfFiles, SIGNAL(triggered()), _signal_mapper, SLOT(map()));
	_signal_mapper->setMapping(netcdfFiles, ImportFileType::NETCDF);
	QAction* petrelFiles = importFiles->addAction("&Petrel Files...");
	connect(petrelFiles, SIGNAL(triggered()), this, SLOT(loadPetrelFiles()));
	QAction* rasterFiles = importFiles->addAction("&Raster Files...");
	connect(rasterFiles, SIGNAL(triggered()), _signal_mapper, SLOT(map()));
	_signal_mapper->setMapping(rasterFiles, ImportFileType::RASTER);
#ifdef OGS_USE_OPENSG
	QAction* rasterPolyFiles = importFiles->addAction("R&aster Files as PolyData...");
	connect(rasterPolyFiles, SIGNAL(triggered()), this, SLOT(map()));
	_signal_mapper->setMapping(rasterPolyFiles, ImportFileType::POLYRASTER);
#endif
#ifdef Shapelib_FOUND
	QAction* shapeFiles = importFiles->addAction("&Shape Files...");
	connect(shapeFiles, SIGNAL(triggered()), _signal_mapper, SLOT(map()));
	_signal_mapper->setMapping(shapeFiles, ImportFileType::SHAPE);
#endif
	QAction* tetgenFiles = importFiles->addAction("&TetGen Files...");
	connect( tetgenFiles, SIGNAL(triggered()), _signal_mapper, SLOT(map()) );
	_signal_mapper->setMapping(tetgenFiles, ImportFileType::TETGEN);

	QAction* vtkFiles = importFiles->addAction("&VTK Files...");
	connect( vtkFiles, SIGNAL(triggered()), _signal_mapper, SLOT(map()) );
	_signal_mapper->setMapping(vtkFiles, ImportFileType::VTK);

	connect(_signal_mapper, SIGNAL(mapped(int)), this, SLOT(open(int)));

	return importFiles;
}

void MainWindow::loadPetrelFiles()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QStringList sfc_file_names = QFileDialog::getOpenFileNames(
	        this, "Select surface data file(s) to import", "", "Petrel files (*)");
	QStringList well_path_file_names = QFileDialog::getOpenFileNames(
	        this, "Select well path data file(s) to import", "", "Petrel files (*)");
	if (sfc_file_names.size() != 0 || well_path_file_names.size() != 0)
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


		QDir dir = QDir(sfc_file_names.at(0));
		settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
	}
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
														settings.value("lastOpenedFileDirectory").toString(),
														"Geosys FEM condition files (*.cnd *.bc *.ic *.st);;All files (* *.*)");
	QDir dir = QDir(fileName);
	settings.setValue("lastOpenedFileDirectory", dir.absolutePath());

	if (!fileName.isEmpty())
		this->loadFEMConditionsFromFile(fileName, geoName);
}

void MainWindow::loadFEMConditionsFromFile(const QString &fileName, std::string geoName)
{
	QFileInfo fi(fileName);
	if (fi.suffix().toLower() == "cnd")
	{
		std::vector<FEMCondition*> conditions;
		std::string schemaName(_fileFinder.getPath("OpenGeoSysCond.xsd"));
		XmlCndInterface xml(&_project, schemaName);
		xml.readFile(conditions, fileName);
		this->addFEMConditions(conditions);
	}
}

void MainWindow::addFEMConditions(const std::vector<FEMCondition*> conditions)
{
	if (!conditions.empty())
	{
		for (size_t i = 0; i < conditions.size(); i++)
		{
			bool condition_ok(true);
			if (conditions[i]->getProcessDistributionType() == FiniteElement::DIRECT)
			{
				if (_meshModels->getMesh(conditions[i]->getAssociatedGeometryName()) != NULL) {
					std::vector<MeshLib::Node*> nodes = _meshModels->getMesh(conditions[i]->getAssociatedGeometryName())->getNodes();
					const size_t nPoints(nodes.size());
					std::vector<GeoLib::Point*> *new_points = new std::vector<GeoLib::Point*>(nPoints);
					for (size_t j = 0; j < nPoints; j++)
						(*new_points)[j] = new GeoLib::Point(nodes[j]->getCoords());
					GeoLib::PointVec pnt_vec("MeshNodes", new_points);
					std::vector<GeoLib::Point*> *cond_points = pnt_vec.getSubset(conditions[i]->getDisNodes());
					std::string geo_name = conditions[i]->getGeoName();
					this->_geoModels->addPointVec(cond_points, geo_name);
					conditions[i]->setGeoName(geo_name); // this might have been changed upon inserting it into geo_objects
				} else {
					OGSError::box("Please load an appropriate geometry first", "Error");
					condition_ok = false;
				}
			}
			if (condition_ok) {
				this->_processModel->addCondition(conditions[i]);
			}
		}
	}
}

void MainWindow::writeFEMConditionsToFile(const QString &geoName, const FEMCondition::CondType type, const QString &fileName)
{
	QFileInfo fi(fileName);
	if (fi.suffix().compare("cnd") == 0 )
	{
		std::string schemaName(_fileFinder.getPath("OpenGeoSysCond.xsd"));
		XmlCndInterface xml(&_project, schemaName);
		xml.setNameForExport(geoName.toStdString());
		xml.setConditionType(type);
		xml.writeToFile(fileName.toStdString());
	}
}

void MainWindow::writeGeometryToFile(QString gliName, QString fileName)
{
	std::string schemaName(_fileFinder.getPath("OpenGeoSysGLI.xsd"));
	XmlGmlInterface xml(&_project, schemaName);
	xml.setNameForExport(gliName.toStdString());
	xml.writeToFile(fileName.toStdString());
}

void MainWindow::writeStationListToFile(QString listName, QString fileName)
{
	std::string schemaName(_fileFinder.getPath("OpenGeoSysSTN.xsd"));
	XmlStnInterface xml(&_project, schemaName);
	xml.setNameForExport(listName.toStdString());
	xml.writeToFile(fileName.toStdString());
}

void MainWindow::exportBoreholesToGMS(std::string listName,
                                      std::string fileName)
{
	const std::vector<GeoLib::Point*>* stations(_geoModels->getStationVec(listName));
	GMSInterface::writeBoreholesToGMS(stations, fileName);
}

void MainWindow::callGMSH(std::vector<std::string> & selectedGeometries,
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
			if (param4 == -1) { // adaptive meshing selected
				GMSHInterface gmsh_io(*(static_cast<GeoLib::GEOObjects*> (_geoModels)), true,
								FileIO::GMSH::AdaptiveMeshDensity, param2, param3, param1,
								selectedGeometries);
				gmsh_io.setPrecision(20);
				gmsh_io.writeToFile(fileName.toStdString());
			} else { // homogeneous meshing selected
				GMSHInterface gmsh_io(*(static_cast<GeoLib::GEOObjects*> (_geoModels)), true,
								FileIO::GMSH::FixedMeshDensity, param4, param3, param1,
								selectedGeometries);
				gmsh_io.setPrecision(20);
				gmsh_io.writeToFile(fileName.toStdString());
			}


			if (system(NULL) != 0) // command processor available
			{
				std::string gmsh_command("gmsh -2 -algo meshadapt ");
				std::string fname (fileName.toStdString());
				gmsh_command += fname;
				size_t pos (fname.rfind ("."));
				if (pos != std::string::npos)
					fname = fname.substr (0, pos);
				gmsh_command += " -o " + fname + ".msh";
				system(gmsh_command.c_str());
				this->loadFile(ImportFileType::GMSH, fileName.left(fileName.length() - 3).append("msh"));
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
		std::cout << "No geometry information selected..." << std::endl;
}

void MainWindow::showConditionWriterDialog()
{
	ConditionWriterDialog dlg(_geoModels);
	connect(&dlg , SIGNAL(saveFEMConditionsRequested(const QString&, const FEMCondition::CondType, const QString&)),
	        this, SLOT(writeFEMConditionsToFile(const QString&, const FEMCondition::CondType, const QString&)));
	dlg.exec();
}

void MainWindow::showDiagramPrefsDialog(QModelIndex &index)
{
	QString listName;
	GeoLib::Station* stn = _geoModels->getStationModel()->stationFromIndex(
	        index, listName);

	if ((stn->type() == GeoLib::Station::STATION) && stn->getSensorData())
	{
		DiagramPrefsDialog* prefs ( new DiagramPrefsDialog(stn) );
		prefs->setAttribute(Qt::WA_DeleteOnClose);
		prefs->show();
	}
	if (stn->type() == GeoLib::Station::BOREHOLE)
		OGSError::box("No time series data available for borehole.");
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
/* TODO6
void MainWindow::showFileConverterDialog()
{
	OGSFileConverter dlg;
	dlg.exec();
}
*/
void MainWindow::showGeoNameDialog(const std::string &geometry_name, const GeoLib::GEOTYPE object_type, size_t id)
{
	std::string old_name = this->_geoModels->getElementNameByID(geometry_name, object_type, id);
	SetNameDialog dlg(geometry_name, GeoLib::convertGeoTypeToString(object_type), id, old_name);
	connect(&dlg, SIGNAL(requestNameChange(const std::string&, const GeoLib::GEOTYPE, size_t, std::string)),
		this->_geoModels, SLOT(addNameForElement(const std::string&, const GeoLib::GEOTYPE, size_t, std::string)));
	dlg.exec();

	static_cast<GeoTreeModel*>(this->geoTabWidget->treeView->model())->setNameForItem(geometry_name, object_type,
		id,	this->_geoModels->getElementNameByID(geometry_name, object_type, id));
}

void MainWindow::showCondSetupDialog(const std::string &geometry_name, const GeoLib::GEOTYPE object_type, size_t id, bool on_points)
{
	std::string geo_name("");
	if (object_type != GeoLib::INVALID)
		geo_name = this->_geoModels->getElementNameByID(geometry_name, object_type, id);
	else
		geo_name = geometry_name; // in this case this is actually the mesh name

	if (geo_name.empty())
	{
		this->showGeoNameDialog(geometry_name, object_type, id);
		geo_name = this->_geoModels->getElementNameByID(geometry_name, object_type, id);
	}
	// Object should now have a name ... if not, cancel the setup process
	if (geo_name.empty())
		OGSError::box("FEM Condition Setup canceled.");
	else
	{
		if (on_points)
			this->_geoModels->addNameForObjectPoints(geometry_name, object_type, geo_name, geometry_name);

		if (object_type != GeoLib::INVALID)
		{
			FEMConditionSetupDialog dlg(geometry_name, object_type, geo_name, this->_geoModels->getGEOObject(geometry_name, object_type, geo_name), on_points);
			connect(&dlg, SIGNAL(createFEMCondition(std::vector<FEMCondition*>)), this, SLOT(addFEMConditions(std::vector<FEMCondition*>)));
			dlg.exec();
		}
		else
		{
			const MeshLib::Mesh* mesh = _project.getMesh(geo_name);
			FEMConditionSetupDialog dlg(geo_name, mesh);
			connect(&dlg, SIGNAL(createFEMCondition(std::vector<FEMCondition*>)), this, SLOT(addFEMConditions(std::vector<FEMCondition*>)));
			dlg.exec();
		}
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
	connect(&dlg, SIGNAL(requestMeshing(std::vector<std::string> &, size_t, double, double, double, bool)),
	        this, SLOT(callGMSH(std::vector<std::string> &, size_t, double, double, double, bool)));
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
	const double dir[3] = {0, 0, 1};
	const MeshLib::Mesh* mesh = this->_project.getMesh("tb_wo_mat");
	_meshModels->addMesh( MeshLib::MshEditor::getMeshSurface(*mesh, dir) );
	/*
	std::vector<GeoLib::PointWithID*> pnts (MeshLib::MshEditor::getSurfaceNodes(*mesh, dir));
	std::vector<GeoLib::Point*> *sfcpnts = new std::vector<GeoLib::Point*>(pnts.size());
	for (unsigned i=0; i<pnts.size(); ++i)
		(*sfcpnts)[i] = pnts[i];
	std::string name("SurfacePoints");
	this->_geoModels->addPointVec(sfcpnts, name);
	*/
/*
	const std::vector<GeoLib::Polyline*> *lines = this->_geoModels->getPolylineVec("WESS Rivers");	MeshLib::CFEMesh* mesh = const_cast<MeshLib::CFEMesh*>(_project.getMesh("Ammer-Homogen100m-Final"));
	std::vector<size_t> nodes;
	mesh->GetNODOnPLY((*lines)[0], nodes);

	std::vector<GeoLib::Point*> *points = new std::vector<GeoLib::Point*>(nodes.size());
	for (size_t i=0; i<nodes.size(); i++)
	{
		(*points)[i] = new GeoLib::Point(mesh->nod_vector[nodes[i]]->getData());
	}
	std::string testpoints_name("testpoints");
	_geoModels->addPointVec(points, testpoints_name);
*/
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

