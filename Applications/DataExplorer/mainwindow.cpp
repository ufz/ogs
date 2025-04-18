/**
 * \file
 * \author Lars Bilke
 * \date   2009-11-04
 * \brief  Implementation of the MainWindow class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (https://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "mainwindow.h"

#include "BaseLib/Logging.h"

// Qt includes
#include <QDate>
#include <QDesktopWidget>
#include <QFileDialog>
#include <QMessageBox>
#include <QObject>
#include <QScreen>
#include <QSettings>
#include <QSignalMapper>
#ifndef NDEBUG
#include <QElapsedTimer>
#endif  // NDEBUG

// VTK includes
#include <vtkImageAlgorithm.h>
#include <vtkImageData.h>
#include <vtkOBJExporter.h>
#include <vtkRenderer.h>
#include <vtkVRMLExporter.h>

#include "Applications/FileIO/FEFLOW/FEFLOWGeoInterface.h"
#include "Applications/FileIO/FEFLOW/FEFLOWMeshInterface.h"
#include "Applications/FileIO/GMSInterface.h"
#include "Applications/FileIO/Gmsh/GMSHInterface.h"
#include "Applications/FileIO/Gmsh/GmshReader.h"
#include "Applications/FileIO/GocadIO/GocadAsciiReader.h"
#include "Applications/FileIO/Legacy/OGSIOVer4.h"
#include "Applications/FileIO/PetrelInterface.h"
#include "Applications/FileIO/TetGenInterface.h"
#include "Applications/FileIO/XmlIO/Qt/XmlPrjInterface.h"
#include "Applications/Utils/OGSFileConverter/OGSFileConverter.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Histogram.h"
#include "GeoLib/DuplicateGeometry.h"
#include "GeoLib/IO/AsciiRasterInterface.h"
#include "GeoLib/IO/XmlIO/Qt/XmlGmlInterface.h"
#include "GeoLib/IO/XmlIO/Qt/XmlStnInterface.h"
#include "GeoLib/Raster.h"
#include "InfoLib/GitInfo.h"
#include "MeshGeoToolsLib/GeoMapper.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/Legacy/MeshIO.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Vtk/VtkMappedMeshSource.h"
#include "MeshToolsLib/MeshQuality/ElementQualityInterface.h"
#include "MeshToolsLib/MeshSurfaceExtraction.h"
#include "MeshToolsLib/convertMeshToGeo.h"

// Dialogs
#include "DataView/AddFaultsToVoxelGridDialog.h"
#include "DataView/CreateStructuredGridDialog.h"
#include "DataView/DataExplorerSettingsDialog.h"
#include "DataView/DiagramView/DiagramPrefsDialog.h"
#include "DataView/GMSHPrefsDialog.h"
#include "DataView/GeoOnMeshMappingDialog.h"
#include "DataView/Layers2GridDialog.h"
#include "DataView/LicenseDialog.h"
#include "DataView/LineEditDialog.h"
#include "DataView/MergeGeometriesDialog.h"
#include "DataView/MeshAnalysisDialog.h"
#include "DataView/MeshElementRemovalDialog.h"
#include "DataView/MeshQualitySelectionDialog.h"
#include "DataView/TranslateDataDialog.h"
#include "DataView/Vtu2GridDialog.h"
#ifdef OGS_USE_NETCDF
#include "VtkVis/NetCdfConfigureDialog.h"
#endif  // OGS_USE_NETCDF
#include "Base/LastSavedFileDirectory.h"
#include "Base/OGSError.h"
#include "Base/RecentFiles.h"
#include "Base/TreeModelIterator.h"
#include "DataView/GeoTreeModel.h"
#include "DataView/SHPImportDialog.h"
#include "DataView/SetNameDialog.h"
#include "DataView/StationTreeModel.h"
#include "VtkVis/VtkAddFilterDialog.h"
#include "VtkVis/VtkBGImageSource.h"
#include "VtkVis/VtkGeoImageSource.h"
#include "VtkVis/VtkRaster.h"
#include "VtkVis/VtkVisPipelineItem.h"

using namespace FileIO;

MainWindow::MainWindow(QWidget* parent /* = 0*/) : QMainWindow(parent)
{
    setupUi(this);

    // Setup various models
    _geo_model = std::make_unique<GEOModels>(_project.getGEOObjects());
    _meshModel = std::make_unique<MeshModel>(_project);
    _elementModel = std::make_unique<ElementTreeModel>();
    _processModel = std::make_unique<ProcessModel>(_project);
    _conditionModel = std::make_unique<FemConditionModel>();

    geoTabWidget->treeView->setModel(_geo_model->getGeoModel());
    stationTabWidget->treeView->setModel(_geo_model->getStationModel());
    meshTabWidget->treeView->setModel(_meshModel.get());
    meshTabWidget->elementView->setModel(_elementModel.get());
    modellingTabWidget->treeView->setModel(_processModel.get());
    modellingTabWidget->conditionView->setModel(_conditionModel.get());

    // vtk visualization pipeline
    _vtkVisPipeline =
        std::make_unique<VtkVisPipeline>(visualizationWidget->renderer());

    // station model connects
    connect(stationTabWidget->treeView, SIGNAL(openStationListFile(int)), this,
            SLOT(open(int)));
    connect(stationTabWidget->treeView,
            SIGNAL(stationListExportRequested(std::string, std::string)), this,
            SLOT(exportBoreholesToGMS(
                std::string, std::string)));  // export Stationlist to GMS
    connect(stationTabWidget->treeView, SIGNAL(stationListRemoved(std::string)),
            _geo_model.get(),
            SLOT(removeStationVec(
                std::string)));  // update model when stations are removed
    connect(stationTabWidget->treeView,
            SIGNAL(stationListSaved(QString, QString)), this,
            SLOT(writeStationListToFile(QString,
                                        QString)));  // save Stationlist to File
    connect(
        _geo_model.get(),
        SIGNAL(stationVectorRemoved(StationTreeModel*, std::string)), this,
        SLOT(updateDataViews()));  // update data view when stations are removed
    connect(stationTabWidget->treeView,
            SIGNAL(requestNameChangeDialog(const std::string&, std::size_t)),
            this, SLOT(showStationNameDialog(const std::string&, std::size_t)));
    connect(stationTabWidget->treeView,
            SIGNAL(geometryMappingRequested(const std::string&)), this,
            SLOT(mapGeometry(const std::string&)));
    connect(stationTabWidget->treeView, SIGNAL(diagramRequested(QModelIndex&)),
            this,
            SLOT(showDiagramPrefsDialog(
                QModelIndex&)));  // connect treeview to diagramview

    // geo model connects
    connect(geoTabWidget->treeView, SIGNAL(openGeometryFile(int)), this,
            SLOT(open(int)));
    connect(geoTabWidget->treeView,
            SIGNAL(listRemoved(std::string, GeoLib::GEOTYPE)), _geo_model.get(),
            SLOT(removeGeometry(std::string, GeoLib::GEOTYPE)));
    connect(geoTabWidget->treeView,
            SIGNAL(geometryMappingRequested(const std::string&)), this,
            SLOT(mapGeometry(const std::string&)));
    connect(
        geoTabWidget->treeView, SIGNAL(saveToFileRequested(QString, QString)),
        this,
        SLOT(writeGeometryToFile(QString, QString)));  // save geometry to file
    connect(geoTabWidget->treeView,
            SIGNAL(requestPointToStationConversion(std::string const&)), this,
            SLOT(convertPointsToStations(std::string const&)));
    connect(
        geoTabWidget->treeView,
        SIGNAL(requestLineEditDialog(const std::string&)), this,
        SLOT(showLineEditDialog(const std::string&)));  // open line edit dialog
    connect(geoTabWidget->treeView,
            SIGNAL(requestNameChangeDialog(const std::string&,
                                           const GeoLib::GEOTYPE, std::size_t)),
            this,
            SLOT(showGeoNameDialog(const std::string&, const GeoLib::GEOTYPE,
                                   std::size_t)));
    connect(_geo_model.get(),
            SIGNAL(geoDataAdded(GeoTreeModel*, std::string, GeoLib::GEOTYPE)),
            this, SLOT(updateDataViews()));
    connect(_geo_model.get(),
            SIGNAL(geoDataRemoved(GeoTreeModel*, std::string, GeoLib::GEOTYPE)),
            this, SLOT(updateDataViews()));
    connect(geoTabWidget->treeView,
            SIGNAL(geoItemSelected(const vtkPolyDataAlgorithm*, int)),
            _vtkVisPipeline.get(),
            SLOT(highlightGeoObject(const vtkPolyDataAlgorithm*, int)));
    connect(geoTabWidget->treeView, SIGNAL(removeGeoItemSelection()),
            _vtkVisPipeline.get(), SLOT(removeHighlightedGeoObject()));
    connect(stationTabWidget->treeView,
            SIGNAL(geoItemSelected(const vtkPolyDataAlgorithm*, int)),
            _vtkVisPipeline.get(),
            SLOT(highlightGeoObject(const vtkPolyDataAlgorithm*, int)));
    connect(stationTabWidget->treeView, SIGNAL(removeGeoItemSelection()),
            _vtkVisPipeline.get(), SLOT(removeHighlightedGeoObject()));

    // Setup connections for mesh models to GUI
    connect(meshTabWidget->treeView, SIGNAL(openMeshFile(int)), this,
            SLOT(open(int)));
    connect(meshTabWidget->treeView,
            SIGNAL(requestMeshRemoval(const QModelIndex&)), _meshModel.get(),
            SLOT(removeMesh(const QModelIndex&)));
    connect(meshTabWidget->treeView,
            SIGNAL(requestMeshRemoval(const QModelIndex&)), _elementModel.get(),
            SLOT(clearView()));
    connect(
        meshTabWidget->treeView,
        SIGNAL(qualityCheckRequested(MeshLib::VtkMappedMeshSource*)),
        this,
        SLOT(showMeshQualitySelectionDialog(MeshLib::VtkMappedMeshSource*)));
    connect(meshTabWidget->treeView,
            SIGNAL(requestMeshToGeometryConversion(const MeshLib::Mesh*)), this,
            SLOT(convertMeshToGeometry(const MeshLib::Mesh*)));
    connect(meshTabWidget->treeView,
            SIGNAL(elementSelected(vtkUnstructuredGridAlgorithm const* const,
                                   unsigned, bool)),
            _vtkVisPipeline.get(),
            SLOT(highlightMeshComponent(
                vtkUnstructuredGridAlgorithm const* const, unsigned, bool)));
    connect(meshTabWidget->treeView, SIGNAL(meshSelected(MeshLib::Mesh const&)),
            this->_elementModel.get(), SLOT(setMesh(MeshLib::Mesh const&)));
    connect(meshTabWidget->treeView, SIGNAL(meshSelected(MeshLib::Mesh const&)),
            meshTabWidget->elementView, SLOT(updateView()));
    connect(
        meshTabWidget->treeView,
        SIGNAL(elementSelected(vtkUnstructuredGridAlgorithm const* const,
                               unsigned, bool)),
        this->_elementModel.get(),
        SLOT(setElement(vtkUnstructuredGridAlgorithm const* const, unsigned)));
    connect(meshTabWidget->treeView,
            SIGNAL(elementSelected(vtkUnstructuredGridAlgorithm const* const,
                                   unsigned, bool)),
            meshTabWidget->elementView, SLOT(updateView()));
    connect(meshTabWidget->treeView,
            SIGNAL(elementSelected(vtkUnstructuredGridAlgorithm const* const,
                                   unsigned, bool)),
            reinterpret_cast<QObject*>(visualizationWidget->interactorStyle()),
            SLOT(removeHighlightActor()));
    connect(meshTabWidget->treeView, SIGNAL(removeSelectedMeshComponent()),
            _vtkVisPipeline.get(), SLOT(removeHighlightedMeshComponent()));
    connect(meshTabWidget->elementView,
            SIGNAL(nodeSelected(vtkUnstructuredGridAlgorithm const* const,
                                unsigned, bool)),
            reinterpret_cast<QObject*>(visualizationWidget->interactorStyle()),
            SLOT(removeHighlightActor()));
    connect(meshTabWidget->elementView,
            SIGNAL(nodeSelected(vtkUnstructuredGridAlgorithm const* const,
                                unsigned, bool)),
            _vtkVisPipeline.get(),
            SLOT(highlightMeshComponent(
                vtkUnstructuredGridAlgorithm const* const, unsigned, bool)));
    connect(meshTabWidget->elementView, SIGNAL(removeSelectedMeshComponent()),
            _vtkVisPipeline.get(), SLOT(removeHighlightedMeshComponent()));

    // Connection for process model to GUI
    connect(modellingTabWidget->treeView,
            SIGNAL(processVarRemoved(QString const&)), _processModel.get(),
            SLOT(removeProcessVariable(QString const&)));
    connect(modellingTabWidget->treeView,
            SIGNAL(conditionRemoved(QString const&, QString const&)),
            _processModel.get(),
            SLOT(removeCondition(QString const&, QString const&)));
    connect(modellingTabWidget->treeView, SIGNAL(clearConditionView()),
            _conditionModel.get(), SLOT(clearView()));
    connect(modellingTabWidget->treeView,
            SIGNAL(processVarSelected(DataHolderLib::FemCondition*)),
            _conditionModel.get(),
            SLOT(setProcessVariable(DataHolderLib::FemCondition*)));
    connect(modellingTabWidget->treeView,
            SIGNAL(conditionSelected(DataHolderLib::FemCondition*)),
            _conditionModel.get(),
            SLOT(setFemCondition(DataHolderLib::FemCondition*)));
    connect(modellingTabWidget->treeView,
            SIGNAL(processVarSelected(DataHolderLib::FemCondition*)),
            modellingTabWidget->conditionView, SLOT(updateView()));
    connect(modellingTabWidget->treeView,
            SIGNAL(conditionSelected(DataHolderLib::FemCondition*)),
            modellingTabWidget->conditionView, SLOT(updateView()));

    // VisPipeline Connects
    connect(_geo_model.get(),
            SIGNAL(geoDataAdded(GeoTreeModel*, std::string, GeoLib::GEOTYPE)),
            _vtkVisPipeline.get(),
            SLOT(addPipelineItem(GeoTreeModel*, std::string, GeoLib::GEOTYPE)));
    connect(
        _geo_model.get(),
        SIGNAL(geoDataRemoved(GeoTreeModel*, std::string, GeoLib::GEOTYPE)),
        _vtkVisPipeline.get(),
        SLOT(removeSourceItem(GeoTreeModel*, std::string, GeoLib::GEOTYPE)));

    connect(_geo_model.get(),
            SIGNAL(stationVectorAdded(StationTreeModel*, std::string)),
            _vtkVisPipeline.get(),
            SLOT(addPipelineItem(StationTreeModel*, std::string)));
    connect(_geo_model.get(),
            SIGNAL(stationVectorRemoved(StationTreeModel*, std::string)),
            _vtkVisPipeline.get(),
            SLOT(removeSourceItem(StationTreeModel*, std::string)));

    connect(_meshModel.get(), SIGNAL(meshAdded(MeshModel*, QModelIndex)),
            _vtkVisPipeline.get(),
            SLOT(addPipelineItem(MeshModel*, QModelIndex)));
    connect(_meshModel.get(), SIGNAL(meshRemoved(MeshModel*, QModelIndex)),
            _vtkVisPipeline.get(),
            SLOT(removeSourceItem(MeshModel*, QModelIndex)));

    connect(_vtkVisPipeline.get(), SIGNAL(vtkVisPipelineChanged()),
            visualizationWidget->vtkWidget, SLOT(update()));
    connect(_vtkVisPipeline.get(), SIGNAL(vtkVisPipelineChanged()),
            vtkVisTabWidget->vtkVisPipelineView, SLOT(expandAll()));
    connect(_vtkVisPipeline.get(), SIGNAL(itemSelected(const QModelIndex&)),
            vtkVisTabWidget->vtkVisPipelineView,
            SLOT(selectItem(const QModelIndex&)));

    vtkVisTabWidget->vtkVisPipelineView->setModel(_vtkVisPipeline.get());
    connect(vtkVisTabWidget->vtkVisPipelineView,
            SIGNAL(requestRemovePipelineItem(QModelIndex)),
            _vtkVisPipeline.get(), SLOT(removePipelineItem(QModelIndex)));
    connect(vtkVisTabWidget->vtkVisPipelineView,
            SIGNAL(requestAddPipelineFilterItem(QModelIndex)), this,
            SLOT(showAddPipelineFilterItemDialog(QModelIndex)));
    connect(vtkVisTabWidget, SIGNAL(requestViewUpdate()), visualizationWidget,
            SLOT(updateView()));

    connect(vtkVisTabWidget->vtkVisPipelineView,
            SIGNAL(actorSelected(vtkProp3D*)),
            reinterpret_cast<QObject*>(visualizationWidget->interactorStyle()),
            SLOT(highlightActor(vtkProp3D*)));
    connect(_vtkVisPipeline.get(), SIGNAL(vtkVisPipelineChanged()),
            visualizationWidget, SLOT(updateView()));

    // Propagates selected vtk object in the pipeline to the pick interactor
    connect(vtkVisTabWidget->vtkVisPipelineView,
            SIGNAL(dataObjectSelected(vtkDataObject*)),
            reinterpret_cast<QObject*>(visualizationWidget->interactorStyle()),
            SLOT(pickableDataObject(vtkDataObject*)));
    connect(reinterpret_cast<QObject*>(visualizationWidget->vtkPickCallback()),
            SIGNAL(actorPicked(vtkProp3D*)),
            vtkVisTabWidget->vtkVisPipelineView, SLOT(selectItem(vtkProp3D*)));
    connect(reinterpret_cast<QObject*>(visualizationWidget->interactorStyle()),
            SIGNAL(elementPicked(vtkUnstructuredGridAlgorithm const* const,
                                 const unsigned)),
            this->_elementModel.get(),
            SLOT(setElement(vtkUnstructuredGridAlgorithm const* const,
                            const unsigned)));
    connect(reinterpret_cast<QObject*>(visualizationWidget->interactorStyle()),
            SIGNAL(elementPicked(vtkUnstructuredGridAlgorithm const* const,
                                 const unsigned)),
            meshTabWidget->elementView, SLOT(updateView()));
    connect(reinterpret_cast<QObject*>(visualizationWidget->interactorStyle()),
            SIGNAL(clearElementView()), this->_elementModel.get(),
            SLOT(clearView()));
    connect(reinterpret_cast<QObject*>(visualizationWidget->interactorStyle()),
            SIGNAL(elementPicked(vtkUnstructuredGridAlgorithm const* const,
                                 const unsigned)),
            this->_vtkVisPipeline.get(),
            SLOT(removeHighlightedMeshComponent()));

    connect(vtkVisTabWidget->vtkVisPipelineView,
            SIGNAL(meshAdded(MeshLib::Mesh*)), _meshModel.get(),
            SLOT(addMesh(MeshLib::Mesh*)));

    // Stack the data dock widgets together
    tabifyDockWidget(geoDock, mshDock);
    tabifyDockWidget(mshDock, modellingDock);
    tabifyDockWidget(modellingDock, stationDock);

    // Restore window geometry
    readSettings();

    // Get info on screens geometry(ies)
    _vtkWidget.reset(visualizationWidget->vtkWidget);
    for (auto const& screen : QGuiApplication::screens())
    {
        _screenGeometries.push_back(screen->availableGeometry());
    }

    // Setup import files menu
    QMenu* import_files_menu = createImportFilesMenu();  // owned by MainWindow
    menu_File->insertMenu(action_Exit, import_files_menu);

    // Setup recent files menu
    RecentFiles* recentFiles =
        new RecentFiles(this, SLOT(openRecentFile()), "recentFileList");
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
            SLOT(showConditionDockWidget(bool)));
    menuWindows->addAction(showModellingDockAction);

    QAction* showVisDockAction = vtkVisDock->toggleViewAction();
    showVisDockAction->setStatusTip(tr("Shows / hides the VTK Pipeline view"));
    connect(showVisDockAction, SIGNAL(triggered(bool)), this,
            SLOT(showVisDockWidget(bool)));
    menuWindows->addAction(showVisDockAction);

    // Presentation mode
    auto* presentationMenu = new QMenu(this);
    presentationMenu->setTitle("Presentation on");
    connect(presentationMenu, SIGNAL(aboutToShow()), this,
            SLOT(createPresentationMenu()));
    menuWindows->insertMenu(showVisDockAction, presentationMenu);

    _visPrefsDialog = std::make_unique<VisPrefsDialog>(*_vtkVisPipeline,
                                                       *visualizationWidget);
}

void MainWindow::closeEvent(QCloseEvent* event)
{
    writeSettings();
    QWidget::closeEvent(event);
}

void MainWindow::showGeoDockWidget(bool show)
{
    if (show)
    {
        geoDock->show();
    }
    else
    {
        geoDock->hide();
    }
}

void MainWindow::showStationDockWidget(bool show)
{
    if (show)
    {
        stationDock->show();
    }
    else
    {
        stationDock->hide();
    }
}

void MainWindow::showMshDockWidget(bool show)
{
    if (show)
    {
        mshDock->show();
    }
    else
    {
        mshDock->hide();
    }
}

void MainWindow::showConditionDockWidget(bool show)
{
    if (show)
    {
        modellingDock->show();
    }
    else
    {
        modellingDock->hide();
    }
}

void MainWindow::showVisDockWidget(bool show)
{
    if (show)
    {
        vtkVisDock->show();
    }
    else
    {
        vtkVisDock->hide();
    }
}

void MainWindow::open(int file_type)
{
    QSettings settings;
    auto t = static_cast<ImportFileType::type>(file_type);
    QString type_str = QString::fromStdString(
        (ImportFileType::convertImportFileTypeToString(t)));
    QString fileName = QFileDialog::getOpenFileName(
        this, "Select " + type_str + " file to import",
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
    auto* action = qobject_cast<QAction*>(sender());
    if (action)
    {
        loadFile(ImportFileType::OGS, action->data().toString());
    }
}

void MainWindow::save()
{
    QString fileName = QFileDialog::getSaveFileName(
        this,
        "Save data as",
        LastSavedFileDirectory::getDir(),
        "GeoSys project (*.prj);;GMSH geometry files (*.geo)");

    if (fileName.isEmpty())
    {
        OGSError::box("No filename specified.");
        return;
    }

    QFileInfo fi(fileName);
    LastSavedFileDirectory::setDir(fileName);

    if (fi.suffix().toLower() == "prj")
    {
        XmlPrjInterface xml(_project);
        BaseLib::IO::writeStringToFile(xml.writeToString(),
                                       fileName.toStdString());
    }
    else if (fi.suffix().toLower() == "geo")
    {
        auto const selected_geometries =
            _project.getGEOObjects().getGeometryNames();

        // values necessary also for the adaptive meshing
        const double point_density = 0;
        const double station_density = point_density;
        const int max_pnts_per_leaf = 0;

        FileIO::GMSH::GMSHInterface gmsh_io(
            _project.getGEOObjects(), true,
            FileIO::GMSH::MeshDensityAlgorithm::FixedMeshDensity, point_density,
            station_density, max_pnts_per_leaf, selected_geometries, false,
            false);
        bool const success = BaseLib::IO::writeStringToFile(
            gmsh_io.writeToString(), fileName.toStdString());

        if (!success)
        {
            OGSError::box(" No geometry available\n to write to geo-file");
        }
    }
}

void MainWindow::loadFile(ImportFileType::type t, const QString& fileName)
{
    QFile file(fileName);
    if (!file.exists())
    {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot read file %1:\n%2.")
                                 .arg(fileName)
                                 .arg(file.errorString()));
        return;
    }

    QApplication::setOverrideCursor(Qt::WaitCursor);
    QFileInfo fi(fileName);
    QSettings settings;
    QDir dir = QDir(fileName);
    std::string base =
        fi.absoluteDir().absoluteFilePath(fi.completeBaseName()).toStdString();

    if (t == ImportFileType::OGS || t == ImportFileType::OGS_GEO ||
        t == ImportFileType::OGS_STN || t == ImportFileType::OGS_MSH)
    {
        if (fi.suffix().toLower() == "gli")
        {
            std::string unique_name;
            std::vector<std::string> errors;
            std::string const gmsh_path =
                settings.value("DataExplorerGmshPath").toString().toStdString();
            if (!FileIO::Legacy::readGLIFileV4(fileName.toStdString(),
                                               _project.getGEOObjects(),
                                               unique_name, errors, gmsh_path))
            {
                for (auto& error : errors)
                {
                    OGSError::box(QString::fromStdString(error));
                }
            }
        }
        else if (fi.suffix().toLower() == "prj")
        {
            XmlPrjInterface xml(_project);
            if (xml.readFile(fileName))
            {
                _meshModel->updateModel();
                _processModel->updateModel();
            }
            else
            {
                OGSError::box(
                    "Failed to load project file.\n Please see console for "
                    "details.");
            }
        }
        else if (fi.suffix().toLower() == "gml")
        {
            GeoLib::IO::XmlGmlInterface xml(_project.getGEOObjects());
            try
            {
                if (!xml.readFile(fileName))
                {
                    OGSError::box(
                        "Failed to load geometry.\n Please see console for "
                        "details.");
                }
            }
            catch (std::runtime_error const& err)
            {
                OGSError::box(err.what(),
                              "Failed to read file `" + fileName + "'");
            }
        }
        // OpenGeoSys observation station files (incl. boreholes)
        else if (fi.suffix().toLower() == "stn")
        {
            GeoLib::IO::XmlStnInterface xml(_project.getGEOObjects());
            if (!xml.readFile(fileName))
            {
                OGSError::box(
                    "Failed to load station data.\n Please see console for "
                    "details.");
            }
        }
        // OpenGeoSys mesh files
        else if (fi.suffix().toLower() == "msh" ||
                 fi.suffix().toLower() == "vtu" ||
                 fi.suffix().toLower() == "vtk")
        {
#ifndef NDEBUG
            QElapsedTimer myTimer;
            myTimer.start();
#endif
            std::unique_ptr<MeshLib::Mesh> mesh(
                MeshLib::IO::readMeshFromFile(fileName.toStdString()));
#ifndef NDEBUG
            INFO("Mesh loading time: {:d} ms.", myTimer.restart());
#endif
            if (mesh)
            {
                _meshModel->addMesh(std::move(mesh));
            }
            else
            {
                OGSError::box("Failed to load mesh file.");
            }
#ifndef NDEBUG
            INFO("Mesh model setup time: {:d} ms.", myTimer.elapsed());
#endif
        }

        settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
        emit fileUsed(fileName);
    }
    else if (t == ImportFileType::FEFLOW)
    {
        if (fi.suffix().toLower() == "fem")  // FEFLOW model files
        {
            FileIO::FEFLOWMeshInterface feflowMeshIO;
            std::unique_ptr<MeshLib::Mesh> mesh(
                feflowMeshIO.readFEFLOWFile(fileName.toStdString()));
            if (mesh)
            {
                _meshModel->addMesh(std::move(mesh));
            }
            else
            {
                OGSError::box("Failed to load a FEFLOW mesh.");
            }
            FileIO::FEFLOWGeoInterface feflowGeoIO;
            feflowGeoIO.readFEFLOWFile(fileName.toStdString(),
                                       _project.getGEOObjects());
        }
        settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
    }
    else if (t == ImportFileType::GMS)
    {
        if (fi.suffix().toLower() == "txt")  // GMS borehole files
        {
            std::vector<GeoLib::Point*> boreholes;
            std::string name = fi.baseName().toStdString();

            if (GMSInterface::readBoreholesFromGMS(boreholes,
                                                   fileName.toStdString()))
            {
                _project.getGEOObjects().addStationVec(std::move(boreholes),
                                                       name);
            }
            else
            {
                OGSError::box("Error reading GMS file.");
            }
        }
        else if (fi.suffix().toLower() == "2dm" ||
                 fi.suffix().toLower() == "3dm")  // GMS mesh files
        {
            std::string name = fileName.toStdString();
            std::unique_ptr<MeshLib::Mesh> mesh(GMSInterface::readMesh(name));
            if (mesh)
            {
                _meshModel->addMesh(std::move(mesh));
            }
            else
            {
                OGSError::box("Failed to load a GMS mesh.");
            }
        }
        settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
    }
    else if (t == ImportFileType::GMSH)
    {
        std::string msh_name(fileName.toStdString());
        if (FileIO::GMSH::isGMSHMeshFile(msh_name))
        {
            std::unique_ptr<MeshLib::Mesh> mesh(FileIO::GMSH::readGMSHMesh(
                msh_name, false /*is_created_with_gmsh2*/));
            if (mesh)
            {
                _meshModel->addMesh(std::move(mesh));
            }
            else
            {
                OGSError::box("Failed to load a GMSH mesh.");
            }
        }
        settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
    }
    else if (t == ImportFileType::GOCAD_TSURF)
    {
        std::string file_name(fileName.toStdString());
        std::vector<std::unique_ptr<MeshLib::Mesh>> meshes;
        if (FileIO::Gocad::GocadAsciiReader::readFile(file_name, meshes))
        {
            for (auto& mesh : meshes)
            {
                if (mesh != nullptr)
                {
                    _meshModel->addMesh(std::move(mesh));
                }
            }
        }
        else
        {
            OGSError::box("Error reading file.");
        }
        settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
    }
#ifdef OGS_USE_NETCDF
    else if (t == ImportFileType::NETCDF)  // CH  01.2012
    {
        NetCdfConfigureDialog dlg(fileName.toStdString());
        dlg.exec();
        if (dlg.getMesh())
        {
            std::unique_ptr<MeshLib::Mesh> mesh(dlg.getMesh());
            mesh->setName(dlg.getName());
            _meshModel->addMesh(std::move(mesh));
        }
        if (dlg.getRaster())
            _vtkVisPipeline->addPipelineItem(dlg.getRaster());

        settings.setValue("lastOpenedRasterFileDirectory", dir.absolutePath());
    }
#endif  // OGS_USE_NETCDF
    else if (t == ImportFileType::RASTER)
    {
        VtkGeoImageSource* geoImage = VtkGeoImageSource::New();
        if (geoImage->readImage(fileName))
        {
            _vtkVisPipeline->addPipelineItem(geoImage);
        }
        else
        {
            geoImage->Delete();
            OGSError::box("Error reading raster.");
        }
        settings.setValue("lastOpenedRasterFileDirectory", dir.absolutePath());
    }
    else if (t == ImportFileType::POLYRASTER)
    {
        QImage raster;
        vtkImageAlgorithm* img = VtkRaster::loadImage(fileName.toStdString());
        VtkBGImageSource* bg = VtkBGImageSource::New();
        double* origin = img->GetOutput()->GetOrigin();
        bg->SetRaster(img, origin[0], origin[1],
                      img->GetOutput()->GetSpacing()[0]);
        bg->SetName(fileName);
        _vtkVisPipeline->addPipelineItem(bg);
        settings.setValue("lastOpenedRasterFileDirectory", dir.absolutePath());
    }
    else if (t == ImportFileType::SHAPE)
    {
        SHPImportDialog dlg(
            fileName.toStdString(), _project.getGEOObjects(),
            settings.value("DataExplorerGmshPath").toString().toStdString());
        dlg.exec();
        QDir dir = QDir(fileName);
        settings.setValue("lastOpenedShapeFileDirectory", dir.absolutePath());
    }
    else if (t == ImportFileType::TETGEN)
    {
        if (fi.suffix().toLower().compare("poly") == 0 ||
            fi.suffix().toLower().compare("smesh") == 0)
        {
            FileIO::TetGenInterface tetgen;
            tetgen.readTetGenGeometry(fileName.toStdString(),
                                      _project.getGEOObjects());
        }
        else
        {
            settings.setValue("lastOpenedTetgenFileDirectory",
                              QFileInfo(fileName).absolutePath());
            QString element_fname(fi.path() + "/" + fi.completeBaseName() +
                                  ".ele");

            if (!fileName.isEmpty())
            {
                FileIO::TetGenInterface tetgen;
                std::unique_ptr<MeshLib::Mesh> mesh(tetgen.readTetGenMesh(
                    fileName.toStdString(), element_fname.toStdString()));
                if (mesh)
                {
                    _meshModel->addMesh(std::move(mesh));
                }
                else
                {
                    OGSError::box("Failed to load a TetGen mesh.");
                }
            }
        }
    }
    else if (t == ImportFileType::VTK)
    {
        _vtkVisPipeline->loadFromFile(fileName);
        settings.setValue("lastOpenedVtkFileDirectory", dir.absolutePath());
    }

    QApplication::restoreOverrideCursor();
    updateDataViews();
}

void MainWindow::updateDataViews()
{
    visualizationWidget->updateViewOnLoad();
    geoTabWidget->treeView->updateView();
    stationTabWidget->treeView->updateView();
    meshTabWidget->treeView->updateView();
}

void MainWindow::readSettings()
{
    QSettings settings;

    restoreGeometry(settings.value("windowGeometry").toByteArray());
    restoreState(settings.value("windowState").toByteArray());
}

void MainWindow::writeSettings()
{
    QSettings settings;

    settings.setValue("windowGeometry", saveGeometry());
    settings.setValue("windowState", saveState());
}

void MainWindow::showLicense()
{
    LicenseDialog dlg;
    dlg.exec();
}

void MainWindow::about()
{
    QString about(
        "<a href='https://www.opengeosys.org'>www.opengeosys.org</a><br /><br "
        "/>");
    about.append(
        QString("Version: %1<br />")
            .arg(QString::fromStdString(GitInfoLib::GitInfo::ogs_version)));

    about.append(
        QString("Git commit: <a "
                "href='https://github.com/ufz/ogs/commit/%1'>%1</a><br />")
            .arg(QString::fromStdString(
                GitInfoLib::GitInfo::git_version_sha1_short)));
    about.append(QString("Built date: %1<br />")
                     .arg(QDate::currentDate().toString(Qt::ISODate)));

    QMessageBox::about(this, "About OpenGeoSys 6", about);
}

QMenu* MainWindow::createImportFilesMenu()
{
    QMenu* importFiles = new QMenu("&Import Files", this);
    importFiles->addAction("&FEFLOW Files...",
                           [this] { open(ImportFileType::FEFLOW); });
    importFiles->addAction("SMS/G&MS Files...",
                           [this] { open(ImportFileType::GMS); });
    importFiles->addAction("&GMSH Files...",
                           [this] { open(ImportFileType::GMSH); });
    importFiles->addAction("&Gocad TSurface...",
                           [this] { open(ImportFileType::GOCAD_TSURF); });
#ifdef OGS_USE_NETCDF
    importFiles->addAction("&NetCDF Files...",
                           [this] { open(ImportFileType::NETCDF); });
#endif  // OGS_USE_NETCDF
    importFiles->addAction("&Petrel Files...", [this] { loadPetrelFiles(); });
    importFiles->addAction("&Raster Files...",
                           [this] { open(ImportFileType::RASTER); });
    importFiles->addAction("&Shape Files...",
                           [this] { open(ImportFileType::SHAPE); });
    importFiles->addAction("&TetGen Files...",
                           [this] { open(ImportFileType::TETGEN); });
    importFiles->addAction("&VTK Files...",
                           [this] { open(ImportFileType::VTK); });
    return importFiles;
}

void MainWindow::loadPetrelFiles()
{
    QSettings settings;
    QStringList sfc_file_names = QFileDialog::getOpenFileNames(
        this, "Select surface data file(s) to import", "", "Petrel files (*)");
    QStringList well_path_file_names = QFileDialog::getOpenFileNames(
        this, "Select well path data file(s) to import", "",
        "Petrel files (*)");
    if (!sfc_file_names.empty() || !well_path_file_names.empty())
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

        PetrelInterface(sfc_files, well_path_files, unique_str,
                        &_project.getGEOObjects());

        QDir dir = QDir(sfc_file_names.at(0));
        settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
    }
}

void MainWindow::showAddPipelineFilterItemDialog(QModelIndex parentIndex)
{
    VtkAddFilterDialog dlg(*_vtkVisPipeline, parentIndex);
    dlg.exec();
}

void MainWindow::writeGeometryToFile(QString gliName, QString fileName)
{
#ifndef NDEBUG
    QFileInfo fi(fileName);
    if (fi.suffix().toLower() == "gli")
    {
        FileIO::Legacy::writeAllDataToGLIFileV4(fileName.toStdString(),
                                                _project.getGEOObjects());
        return;
    }
#endif
    GeoLib::IO::XmlGmlInterface xml(_project.getGEOObjects());
    xml.export_name = gliName.toStdString();
    BaseLib::IO::writeStringToFile(xml.writeToString(), fileName.toStdString());
}

void MainWindow::writeStationListToFile(QString listName, QString fileName)
{
    GeoLib::IO::XmlStnInterface xml(_project.getGEOObjects());
    xml.export_name = listName.toStdString();
    BaseLib::IO::writeStringToFile(xml.writeToString(), fileName.toStdString());
}

void MainWindow::mapGeometry(const std::string& geo_name)
{
    GeoOnMeshMappingDialog dlg(this->_project.getMeshObjects());
    if (dlg.exec() != QDialog::Accepted)
    {
        return;
    }

    int choice(dlg.getDataSetChoice());

    QString file_name("");
    if (choice < 2)  // load something from a file
    {
        QString file_type[2] = {"OpenGeoSys mesh files (*.vtu *.msh)",
                                "Raster files(*.asc *.grd *.xyz)"};
        QSettings settings;
        file_name = QFileDialog::getOpenFileName(
            this,
            "Select file for mapping",
            settings.value("lastOpenedFileDirectory").toString(),
            file_type[choice]);
        if (file_name.isEmpty())
        {
            return;
        }
        QDir dir = QDir(file_name);
        settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
    }

    MeshGeoToolsLib::GeoMapper geo_mapper(_project.getGEOObjects(), geo_name);
    QFileInfo fi(file_name);
    if (choice == 1)  // load raster from file
    {
        if (fi.suffix().toLower() == "asc" || fi.suffix().toLower() == "grd")
        {
            std::unique_ptr<GeoLib::Raster> raster(
                FileIO::AsciiRasterInterface::getRasterFromASCFile(
                    file_name.toStdString()));
            if (raster)
            {
                geo_mapper.mapOnDEM(std::move(raster));
            }
            else
            {
                OGSError::box("Error reading raster file.");
            }
            _geo_model->updateGeometry(geo_name);
        }
        else
        {
            OGSError::box("The selected file is no supported raster file.");
        }
        return;
    }

    MeshLib::Mesh* mesh(nullptr);
    if (choice == 0)  // load mesh from file
    {
        if (fi.suffix().toLower() == "vtu" || fi.suffix().toLower() == "vtk" ||
            fi.suffix().toLower() == "msh")
        {
            mesh = MeshLib::IO::readMeshFromFile(file_name.toStdString());
        }
        else
        {
            OGSError::box("The selected file is no supported mesh file.");
            return;
        }
    }
    else
    {  // use mesh from ProjectData
        mesh = _project.getMeshObjects()[choice - 2].get();
    }

    std::string new_geo_name = dlg.getNewGeoName();

    if (new_geo_name.empty())
    {
        geo_mapper.mapOnMesh(mesh);
        _geo_model->updateGeometry(geo_name);
    }
    else
    {
        GeoLib::DuplicateGeometry dup(_project.getGEOObjects(), geo_name,
                                      new_geo_name);
        new_geo_name = dup.getFinalizedOutputName();
        MeshGeoToolsLib::GeoMapper mapper(_project.getGEOObjects(),
                                          new_geo_name);
        mapper.advancedMapOnMesh(*mesh);
        _geo_model->updateGeometry(new_geo_name);
    }
    if (choice == 0)
    {
        delete mesh;
    }
}

void MainWindow::convertMeshToGeometry(const MeshLib::Mesh* mesh)
{
    MeshToolsLib::convertMeshToGeo(*mesh, _project.getGEOObjects());
}

void MainWindow::exportBoreholesToGMS(std::string listName,
                                      std::string fileName)
{
    const std::vector<GeoLib::Point*>* stations(
        _project.getGEOObjects().getStationVec(listName));
    GMSInterface::writeBoreholesToGMS(stations, fileName);
}

void MainWindow::callGMSH(std::vector<std::string>& selectedGeometries,
                          unsigned param1, double param2, double param3,
                          double param4, bool delete_geo_file)
{
    if (!selectedGeometries.empty())
    {
        INFO("Start meshing ...");

        QString fileName("");
        QString dir_str = this->getLastUsedDir();

        if (!delete_geo_file)
        {
            fileName = QFileDialog::getSaveFileName(
                this, "Save GMSH-file as",
                LastSavedFileDirectory::getDir() + "tmp_gmsh.geo",
                "GMSH geometry files (*.geo)");
        }
        else
        {
            fileName = "tmp_gmsh.geo";
        }

        if (!fileName.isEmpty())
        {
            try
            {
                if (param4 == -1)
                {  // adaptive meshing selected
                    FileIO::GMSH::GMSHInterface gmsh_io(
                        _project.getGEOObjects(), true,
                        FileIO::GMSH::MeshDensityAlgorithm::AdaptiveMeshDensity,
                        param2, param3, param1, selectedGeometries, false,
                        false);
                    BaseLib::IO::writeStringToFile(gmsh_io.writeToString(),
                                                   fileName.toStdString());
                }
                else
                {  // homogeneous meshing selected
                    FileIO::GMSH::GMSHInterface gmsh_io(
                        _project.getGEOObjects(), true,
                        FileIO::GMSH::MeshDensityAlgorithm::FixedMeshDensity,
                        param4, param3, param1, selectedGeometries, false,
                        false);
                    BaseLib::IO::writeStringToFile(gmsh_io.writeToString(),
                                                   fileName.toStdString());
                }

                if (system(nullptr) != 0)  // command processor available
                {
                    QSettings settings;
                    std::string gmsh_path =
                        settings.value("DataExplorerGmshPath")
                            .toString()
                            .toStdString();

                    if (!gmsh_path.empty())
                    {
                        std::string fname(fileName.toStdString());
                        std::string gmsh_command =
                            "\"" + gmsh_path + "\" -2 -algo meshadapt " + fname;
                        std::size_t pos(fname.rfind("."));
                        if (pos != std::string::npos)
                        {
                            fname = fname.substr(0, pos);
                        }
                        gmsh_command += " -o " + fname + ".msh";
                        // Newer gmsh versions write a newer file format for
                        // meshes per default. At the moment we can't read this
                        // new format. This is a switch for gmsh to write the
                        // 'old' file format.
                        gmsh_command += " -format msh22";
                        auto const return_value =
                            std::system(gmsh_command.c_str());
                        if (return_value != 0)
                        {
                            QString const message =
                                "Execution of gmsh command returned non-zero "
                                "status, " +
                                QString::number(return_value);
                            OGSError::box(message, "Error");
                        }
                        else
                        {
                            this->loadFile(ImportFileType::GMSH,
                                           fileName.left(fileName.length() - 3)
                                               .append("msh"));
                        }
                    }
                    else
                    {
                        OGSError::box("Location of GMSH not specified.",
                                      "Error");
                    }
                }
                else
                {
                    OGSError::box(
                        "Error executing command gmsh - no command processor "
                        "available",
                        "Error");
                }
            }
            catch (std::runtime_error& error)
            {
                OGSError::box(QString(error.what()) +
                                  QString("\n Please cleanup the input data."),
                              "ERROR");
            }
            if (delete_geo_file)
            {
                BaseLib::removeFile(fileName.toStdString());
            }
        }
    }
    else
        INFO("No geometry information selected.");
    QApplication::restoreOverrideCursor();
}

void MainWindow::showFileConverter()
{
    QSettings settings;
    auto* dlg = new OGSFileConverter(
        settings.value("DataExplorerGmshPath").toString().toStdString(), this);
    dlg->setAttribute(Qt::WA_DeleteOnClose);
    dlg->show();
    dlg->raise();
}

void MainWindow::showDiagramPrefsDialog(QModelIndex& index)
{
    QString listName;
    GeoLib::Station* stn =
        _geo_model->getStationModel()->stationFromIndex(index, listName);

    if (dynamic_cast<GeoLib::StationBorehole*>(stn))
    {
        OGSError::box("No time series data available for borehole.");
    }
    else if (dynamic_cast<GeoLib::Station*>(stn) && stn->getSensorData())
    {
        auto* prefs(new DiagramPrefsDialog(stn));
        prefs->setAttribute(Qt::WA_DeleteOnClose);
        prefs->show();
    }
}

void MainWindow::showDiagramPrefsDialog()
{
    QSettings settings;
    QString fileName = QFileDialog::getOpenFileName(
        this, "Select data file to open",
        settings.value("lastOpenedFileDirectory").toString(),
        "Text files (*.txt);;All files (* *.*)");
    if (!fileName.isEmpty())
    {
        QDir dir = QDir(fileName);
        settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
        auto* prefs = new DiagramPrefsDialog(fileName);
        prefs->setAttribute(Qt::WA_DeleteOnClose);
        prefs->show();
    }
}

void MainWindow::showGeoNameDialog(const std::string& geometry_name,
                                   const GeoLib::GEOTYPE object_type,
                                   std::size_t id)
{
    std::string old_name = _project.getGEOObjects().getElementNameByID(
        geometry_name, object_type, id);
    SetNameDialog dlg(GeoLib::convertGeoTypeToString(object_type), id,
                      old_name);
    if (dlg.exec() != QDialog::Accepted)
    {
        return;
    }

    _geo_model->addNameForElement(geometry_name, object_type, id,
                                  dlg.getNewName());
    static_cast<GeoTreeModel*>(this->geoTabWidget->treeView->model())
        ->setNameForItem(geometry_name, object_type, id,
                         _project.getGEOObjects().getElementNameByID(
                             geometry_name, object_type, id));
}

void MainWindow::showStationNameDialog(const std::string& stn_vec_name,
                                       std::size_t id)
{
    std::vector<GeoLib::Point*> const* stations =
        _project.getGEOObjects().getStationVec(stn_vec_name);
    auto* const stn = static_cast<GeoLib::Station*>((*stations)[id]);
    SetNameDialog dlg("Station", id, stn->getName());
    if (dlg.exec() != QDialog::Accepted)
    {
        return;
    }

    stn->setName(dlg.getNewName());
    static_cast<StationTreeModel*>(this->stationTabWidget->treeView->model())
        ->setNameForItem(stn_vec_name, id, stn->getName());
}

void MainWindow::showCreateStructuredGridDialog()
{
    CreateStructuredGridDialog dlg;
    connect(&dlg, SIGNAL(meshAdded(MeshLib::Mesh*)), _meshModel.get(),
            SLOT(addMesh(MeshLib::Mesh*)));
    dlg.exec();
}

void MainWindow::showMeshElementRemovalDialog()
{
    MeshElementRemovalDialog dlg(this->_project);
    connect(&dlg, SIGNAL(meshAdded(MeshLib::Mesh*)), _meshModel.get(),
            SLOT(addMesh(MeshLib::Mesh*)));
    dlg.exec();
}

void MainWindow::showMeshAnalysisDialog()
{
    auto* dlg = new MeshAnalysisDialog(this->_project.getMeshObjects());
    dlg->exec();
}

void MainWindow::showTranslateDataDialog()
{
    auto dlg = TranslateDataDialog(_meshModel.get(), _geo_model.get());
    dlg.exec();
}

void MainWindow::showLayers2GridDialog()
{
    if (_meshModel == nullptr)
    {
        OGSError::box("The given mesh model does not exist.");
    }

    auto dlg = Layers2GridDialog(*_meshModel);
    dlg.exec();
}

void MainWindow::showAddFaultsToVoxelGridDialog()
{
    if (_meshModel == nullptr)
    {
        OGSError::box("The given mesh model does not exist.");
    }
    AddFaultsToVoxelGridDialog(*_meshModel).exec();
}

void MainWindow::showVtu2GridDialog()
{
    if (_meshModel == nullptr)
    {
        OGSError::box("The given mesh model does not exist.");
    }
    auto dlg = Vtu2GridDialog(*_meshModel);
    dlg.exec();
}

void MainWindow::convertPointsToStations(std::string const& geo_name)
{
    std::string stn_name = geo_name + " Stations";
    int ret = GeoLib::geoPointsToStations(_project.getGEOObjects(), geo_name,
                                          stn_name);
    if (ret == 1)
    {
        OGSError::box("No points found to convert.");
    }
}

void MainWindow::showLineEditDialog(const std::string& geoName)
{
    LineEditDialog lineEdit(
        *(_project.getGEOObjects().getPolylineVecObj(geoName)));
    connect(
        &lineEdit,
        SIGNAL(connectPolylines(const std::string&, std::vector<std::size_t>,
                                double, std::string, bool, bool)),
        _geo_model.get(),
        SLOT(connectPolylineSegments(const std::string&,
                                     std::vector<std::size_t>, double,
                                     std::string, bool, bool)));
    lineEdit.exec();
}

void MainWindow::showGMSHPrefsDialog()
{
    GMSHPrefsDialog dlg(_project.getGEOObjects());
    connect(&dlg,
            SIGNAL(requestMeshing(std::vector<std::string>&, unsigned, double,
                                  double, double, bool)),
            this,
            SLOT(callGMSH(std::vector<std::string>&, unsigned, double, double,
                          double, bool)));
    dlg.exec();
}

void MainWindow::showMergeGeometriesDialog()
{
    MergeGeometriesDialog dlg(_project.getGEOObjects());
    if (dlg.exec() != QDialog::Accepted)
    {
        return;
    }
    std::string name(dlg.getGeometryName());
    if (_project.getGEOObjects().mergeGeometries(dlg.getSelectedGeometries(),
                                                 name) < 0)
    {
        OGSError::box("Points are missing for\n at least one geometry.");
    }
}

void MainWindow::showMeshQualitySelectionDialog(
    MeshLib::VtkMappedMeshSource* mshSource)
{
    if (mshSource == nullptr)
    {
        return;
    }

    MeshQualitySelectionDialog dlg;
    if (dlg.exec() != QDialog::Accepted)
    {
        return;
    }
    MeshLib::MeshQualityType const type(dlg.getSelectedMetric());
    MeshToolsLib::ElementQualityInterface quality_interface(
        *mshSource->GetMesh(), type);
    _vtkVisPipeline->showMeshElementQuality(
        mshSource, type, quality_interface.getQualityVector());

    if (dlg.getHistogram())
    {
        quality_interface.writeHistogram(dlg.getHistogramPath());
    }
}

void MainWindow::showVisalizationPrefsDialog()
{
    _visPrefsDialog->show();
}

void MainWindow::showDataExplorerSettingsDialog()
{
    DataExplorerSettingsDialog dlg;
    dlg.exec();
}

void MainWindow::FEMTestStart() {}

void MainWindow::ShowWindow()
{
    this->show();
}

void MainWindow::HideWindow()
{
    this->hide();
}

void MainWindow::loadFileOnStartUp(const QString& fileName)
{
    QString ext = QFileInfo(fileName).suffix();
    if (ext == "msh" || ext == "vtu" || ext == "gli" || ext == "gml")
    {
        this->loadFile(ImportFileType::OGS, fileName);
    }
}

void MainWindow::on_actionExportVTK_triggered(bool checked /*= false*/)
{
    Q_UNUSED(checked)
    QSettings settings;
    int count = 0;
    QString const filename = QFileDialog::getSaveFileName(
        this,
        "Export object to vtk-files",
        settings.value("lastExportedFileDirectory").toString(),
        "VTK files (*.vtp *.vtu)");
    if (!filename.isEmpty())
    {
        QDir const dir = QDir(filename);
        settings.setValue("lastExportedFileDirectory", dir.absolutePath());

        std::string const basename =
            QFileInfo(filename).path().toStdString() + "/" +
            QFileInfo(filename).baseName().toStdString();
        TreeModelIterator it(_vtkVisPipeline.get());
        ++it;
        while (*it)
        {
            std::string const name = basename + std::to_string(++count) + "-" +
                                     (*it)->data(0).toString().toStdString();
            static_cast<VtkVisPipelineItem*>(*it)->writeToFile(name);
            ++it;
        }
    }
}

void MainWindow::on_actionExportVRML2_triggered(bool checked /*= false*/)
{
    Q_UNUSED(checked)
    QSettings settings;
    QString fileName = QFileDialog::getSaveFileName(
        this, "Save scene to VRML file",
        settings.value("lastExportedFileDirectory").toString(),
        "VRML files (*.wrl);;");
    if (!fileName.isEmpty())
    {
        QDir dir = QDir(fileName);
        settings.setValue("lastExportedFileDirectory", dir.absolutePath());

        vtkVRMLExporter* exporter = vtkVRMLExporter::New();
        exporter->SetFileName(fileName.toStdString().c_str());
        exporter->SetRenderWindow(
            visualizationWidget->vtkWidget->renderWindow());
        exporter->Write();
        exporter->Delete();
    }
}

void MainWindow::on_actionExportObj_triggered(bool checked /*= false*/)
{
    Q_UNUSED(checked)
    QSettings settings;
    QString fileName = QFileDialog::getSaveFileName(
        this, "Save scene to Wavefront OBJ files",
        settings.value("lastExportedFileDirectory").toString(), ";;");
    if (!fileName.isEmpty())
    {
        QDir dir = QDir(fileName);
        settings.setValue("lastExportedFileDirectory", dir.absolutePath());

        vtkOBJExporter* exporter = vtkOBJExporter::New();
        exporter->SetFilePrefix(fileName.toStdString().c_str());
        exporter->SetRenderWindow(
            visualizationWidget->vtkWidget->renderWindow());
        exporter->Write();
        exporter->Delete();
    }
}

void MainWindow::createPresentationMenu()
{
    auto* menu = static_cast<QMenu*>(QObject::sender());
    menu->clear();
    if (!_vtkWidget->parent())
    {
        QAction* action = new QAction("Quit presentation mode", menu);
        connect(action, SIGNAL(triggered()), this,
                SLOT(quitPresentationMode()));
        action->setShortcutContext(Qt::WidgetShortcut);
        action->setShortcut(QKeySequence(Qt::Key_Escape));
        menu->addAction(action);
    }
    else
    {
        int count = 0;
        const int currentScreen =
            QApplication::desktop()->screenNumber(visualizationWidget);
        foreach (QRect screenGeo, _screenGeometries)
        {
            Q_UNUSED(screenGeo);
            QAction* action =
                new QAction(QString("On screen %1").arg(count), menu);
            connect(action, SIGNAL(triggered()), this,
                    SLOT(startPresentationMode()));
            if (count == currentScreen)
            {
                action->setEnabled(false);
            }
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
    QString actionText = static_cast<QAction*>(QObject::sender())->text();
    int screen = actionText.split(" ").back().toInt();

    // Move the widget to the screen and maximize it
    // Real fullscreen hides the menu
    _vtkWidget->setParent(nullptr, Qt::Window);
    _vtkWidget->move(
        QPoint(_screenGeometries[screen].x(), _screenGeometries[screen].y()));
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
    visualizationWidget->layout()->addWidget(_vtkWidget.get());

    QMainWindow::centralWidget()->show();

    // Restore the previously saved QMainWindow state
    this->restoreState(_windowState);
}

QString MainWindow::getLastUsedDir()
{
    QSettings settings;
    QString fileName("");
    QStringList files = settings.value("recentFileList").toStringList();
    if (!files.empty())
    {
        return QFileInfo(files[0]).absolutePath();
    }

    return QDir::homePath();
}
