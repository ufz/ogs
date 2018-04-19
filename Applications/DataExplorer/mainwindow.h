/**
 * \file
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "Applications/DataHolderLib/Project.h"

#include "ImportFileTypes.h"
#include "ui_mainwindow.h"

#include "GEOModels.h"
#include "MshModel.h"
#include "ProcessModel.h"
#include "ElementTreeModel.h"
#include "FemConditionModel.h"
#include "VisPrefsDialog.h"
#include "VtkVisPipeline.h"

class TreeModel;

namespace MeshLib
{
    class VtkMappedMeshSource;
}

class QSignalMapper;

/**
 * Main program window for the graphical user interface of OpenGeoSys.
 */
class MainWindow : public QMainWindow, public Ui_MainWindowClass
{
    Q_OBJECT

public:
    MainWindow(QWidget* parent = nullptr);

    void ShowWindow();
    void HideWindow();
    void loadFileOnStartUp(const QString &fileName);

protected:
    void closeEvent(QCloseEvent* event) override;

protected slots:
    void showGeoDockWidget( bool show );
    void showStationDockWidget( bool show );
    void showMshDockWidget( bool show );
    void showConditionDockWidget( bool show );
    void showVisDockWidget( bool show );

    /// Function calls for opening files.
    void open(int i = 0);
    /// Function calls for saving files.
    void save();
    /// Function calls for generating GMSH files from the GUI
    void callGMSH(std::vector<std::string> & selectedGeometries,
                  unsigned param1,
                  double   param2,
                  double   param3,
                  double   param4,
                  bool     delete_geo_file);
    /// Function calls for GMS export.
    void exportBoreholesToGMS(std::string listName, std::string fileName);
    /// Testing functionality for connection to FEM lib
    void FEMTestStart();
    void loadPetrelFiles();
    void mapGeometry(const std::string &geo_name);
    void convertMeshToGeometry(const MeshLib::Mesh* mesh);
    void openRecentFile();
    void about();
    void showAddPipelineFilterItemDialog(QModelIndex parentIndex);
    void showDataExplorerSettingsDialog();
    /// Allows setting the name for a geometric object
    void showGeoNameDialog(const std::string &geometry_name, const GeoLib::GEOTYPE object_type, std::size_t id);
    /// Allows setting the name for a station
    void showStationNameDialog(const std::string& stn_vec_name, std::size_t id);
    /// Creates a structured grid with user-specified parameters.
    void showCreateStructuredGridDialog();
    /// Removal of mesh elements based on a number of criteria.
    void showMeshElementRemovalDialog();
    /// Calls the diagram prefs dialog from the Tools menu.
    void showDiagramPrefsDialog();
    /// Calls the diagram prefs dialog from the station list (i.e. for a specific station).
    void showDiagramPrefsDialog(QModelIndex &index);
    /// Calls the OGSFileConverter as an external tool
    void showFileConverter();
    //TODO6 void showFileConverterDialog();
    void showLicense();
    void showLineEditDialog(const std::string &geoName);
    void showGMSHPrefsDialog();
    void showMergeGeometriesDialog();
    void showMeshAnalysisDialog();
    void showMeshQualitySelectionDialog(MeshLib::VtkMappedMeshSource* mshSource);
    void showVisalizationPrefsDialog();
    void updateDataViews();
    void writeGeometryToFile(QString listName, QString fileName);
    void writeStationListToFile(QString listName, QString fileName);

    void on_actionExportVTK_triggered(bool checked = false);
    void on_actionExportVRML2_triggered(bool checked = false);
    void on_actionExportObj_triggered(bool checked = false);

    void createPresentationMenu();
    void startPresentationMode();
    void quitPresentationMode();

private:
    QMenu* createImportFilesMenu();
    void loadFile(ImportFileType::type t, const QString &fileName);
    void loadFEMConditionsFromFile(const QString &fileName, std::string geoName = "");
    void readSettings();
    void writeSettings();
    QString getLastUsedDir();

    DataHolderLib::Project _project;
    std::unique_ptr<MshModel> _meshModel;
    std::unique_ptr<ElementTreeModel> _elementModel;
    std::unique_ptr<ProcessModel> _processModel;
    std::unique_ptr<FemConditionModel> _conditionModel;
    std::unique_ptr<VtkVisPipeline> _vtkVisPipeline;
    QList<QRect> _screenGeometries;
    std::unique_ptr<QWidget> _vtkWidget;
    QByteArray _windowState;

    std::unique_ptr<VisPrefsDialog> _visPrefsDialog;

    std::unique_ptr<GEOModels> _geo_model;

signals:
    void fileUsed( QString filename );
    void fileOpenRequested( int );
};
