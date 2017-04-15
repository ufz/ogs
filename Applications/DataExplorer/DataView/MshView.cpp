/**
 * \file
 * \author Lars Bilke
 * \date   2009-09-24
 * \brief  Implementation of the MshView class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MshView.h"

#include <memory>

#include <QHeaderView>
#include <QContextMenuEvent>
#include <QFileDialog>
#include <QMenu>
#include <QObject>
#include <QSettings>

#include "SHPInterface.h"
#include "TetGenInterface.h"
#include "Applications/FileIO/AsciiRasterInterface.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/MeshSurfaceExtraction.h"
#include "MeshLib/MeshEditing/AddLayerToMesh.h"

#include "OGSError.h"
#include "MeshMapping2DDialog.h"
#include "MeshLayerEditDialog.h"
#include "MeshValueEditDialog.h"
#include "SurfaceExtractionDialog.h"
#include "AddLayerToMeshDialog.h"
#include "MshItem.h"
#include "MshModel.h"

#include "ImportFileTypes.h"
#include "LastSavedFileDirectory.h"
#include "SaveMeshDialog.h"


MshView::MshView( QWidget* parent /*= 0*/ )
    : QTreeView(parent)
{
    setUniformRowHeights(true);
    //resizeColumnsToContents();
    //resizeRowsToContents();
}

MshView::~MshView() = default;

void MshView::updateView()
{
    setAlternatingRowColors(true);
    setColumnWidth(0,125);
    std::size_t nColumns = (this->model() != nullptr) ? this->model()->columnCount() : 0;
    for (std::size_t i = 1; i < nColumns; i++)
        resizeColumnToContents(i);
}

void MshView::selectionChanged( const QItemSelection &selected, const QItemSelection &deselected )
{
    Q_UNUSED(deselected);
    if (!selected.isEmpty())
    {
        emit removeSelectedMeshComponent();
        const QModelIndex idx = *(selected.indexes().begin());
        const TreeItem* tree_item = static_cast<TreeModel*>(this->model())->getItem(idx);

        const MshItem* list_item = dynamic_cast<const MshItem*>(tree_item);
        if (list_item)
        {
            emit enableSaveButton(true);
            emit enableRemoveButton(true);
            emit meshSelected(*list_item->getMesh());
        }
        else
        {
            emit enableSaveButton(false);
            emit enableRemoveButton(false);
            emit elementSelected(dynamic_cast<const MshItem*>(tree_item->parentItem())->vtkSource(), static_cast<unsigned>(tree_item->row()), true);
        }
    }
}

void MshView::addMesh()
{
    emit openMeshFile(ImportFileType::OGS_MSH);
}

void MshView::removeMesh()
{
    QModelIndex index (this->selectionModel()->currentIndex());
    if (!index.isValid())
        OGSError::box("No mesh selected.");
    else
    {
        emit requestMeshRemoval(index);
        emit enableSaveButton(false);
        emit enableRemoveButton(false);
    }
}

void MshView::contextMenuEvent( QContextMenuEvent* event )
{
    QModelIndex const& index = this->selectionModel()->currentIndex();
    MshItem const*const item = dynamic_cast<MshItem*>(static_cast<TreeItem*>(index.internalPointer()));

    if (item == nullptr)
        return;

    unsigned const mesh_dim (item->getMesh()->getDimension());

    std::vector<MeshAction> actions;
    actions.push_back({new QAction("Map mesh...", this), 1, 2});
    connect(actions.back().action, SIGNAL(triggered()), this, SLOT(openMap2dMeshDialog()));
    actions.push_back({new QAction("Edit mesh...", this), 2, 3});
    connect(actions.back().action, SIGNAL(triggered()), this, SLOT(openMeshEditDialog()));
    actions.push_back({new QAction("Add layer...", this), 1, 3});
    connect(actions.back().action, SIGNAL(triggered()), this, SLOT(openAddLayerDialog()));
    actions.push_back({new QAction("Edit material groups...", this), 1, 3});
    connect(actions.back().action, SIGNAL(triggered()), this, SLOT(openValuesEditDialog()));
    actions.push_back({new QAction("Extract surface...", this), 3, 3});
    connect(actions.back().action, SIGNAL(triggered()), this, SLOT(extractSurfaceMesh()));
    actions.push_back({new QAction("Calculate element quality...", this), 2, 3});
    connect(actions.back().action, SIGNAL(triggered()), this, SLOT(checkMeshQuality()));
    actions.push_back({new QAction("Convert to geometry", this), 1, 2});
    connect(actions.back().action, SIGNAL(triggered()), this, SLOT(convertMeshToGeometry()));
    actions.push_back({new QAction("Export to Shapefile...", this), 2, 2});
    connect(actions.back().action, SIGNAL(triggered()), this, SLOT(exportToShapefile()));
    actions.push_back({new QAction("Export to TetGen...", this), 3, 3});
    connect(actions.back().action, SIGNAL(triggered()), this, SLOT(exportToTetGen()));

    QMenu menu(this);
    for (MeshAction a : actions)
    {
        if (mesh_dim >= a.min_dim && mesh_dim <= a.max_dim)
            menu.addAction(a.action);
    }
    menu.exec(event->globalPos());
}

void MshView::openMap2dMeshDialog()
{
    MshModel const*const model = static_cast<MshModel*>(this->model());
    QModelIndex const index = this->selectionModel()->currentIndex();
    MeshLib::Mesh const*const mesh = model->getMesh(index);
    if (mesh == nullptr)
        return;

    MeshMapping2DDialog dlg;
    if (dlg.exec() != QDialog::Accepted)
        return;

    std::unique_ptr<MeshLib::Mesh> result(new MeshLib::Mesh(*mesh));
    result->setName(dlg.getNewMeshName());
    if (dlg.useRasterMapping())
    {
        std::unique_ptr<GeoLib::Raster> raster{
            FileIO::AsciiRasterInterface::readRaster(dlg.getRasterPath())};
        if (!raster)
        {
            OGSError::box(QString::fromStdString(
                "Error mapping mesh. Could not read raster file " +
                dlg.getRasterPath()));
            return;
        }
        if (!MeshLib::MeshLayerMapper::layerMapping(*result, *raster,
                                                    dlg.getNoDataReplacement()))
        {
            OGSError::box("Error mapping mesh.");
            return;
        }
    }
    else
        MeshLib::MeshLayerMapper::mapToStaticValue(*result, dlg.getStaticValue());
    static_cast<MshModel*>(this->model())->addMesh(std::move(result));

}

void MshView::openMeshEditDialog()
{
    MshModel const*const model = static_cast<MshModel*>(this->model());
    QModelIndex const index = this->selectionModel()->currentIndex();
    MeshLib::Mesh const*const mesh = model->getMesh(index);

    MeshLayerEditDialog meshLayerEdit(mesh);
    connect(&meshLayerEdit, SIGNAL(mshEditFinished(MeshLib::Mesh*)),
            model, SLOT(addMesh(MeshLib::Mesh*)));
    meshLayerEdit.exec();
}

void MshView::openValuesEditDialog()
{
    MshModel const*const model = static_cast<MshModel*>(this->model());
    QModelIndex const index = this->selectionModel()->currentIndex();
    MeshLib::Mesh* mesh = const_cast<MeshLib::Mesh*>(model->getMesh(index));

    MeshValueEditDialog valueEdit(mesh);
    connect(&valueEdit, SIGNAL(valueEditFinished(MeshLib::Mesh*)),
            model, SLOT(updateMesh(MeshLib::Mesh*)));
    valueEdit.exec();
}

void MshView::openAddLayerDialog()
{
    QModelIndex const index = this->selectionModel()->currentIndex();
    if (!index.isValid())
        return;

    MeshLib::Mesh const*const mesh = static_cast<MshModel*>(this->model())->getMesh(index);
    if (mesh == nullptr)
        return;

    AddLayerToMeshDialog dlg;
    if (dlg.exec() != QDialog::Accepted)
        return;

    double const thickness (dlg.getThickness());
    std::unique_ptr<MeshLib::Mesh> result(MeshLib::addLayerToMesh(
        *mesh, thickness, dlg.getName(), dlg.isTopLayer()));

    if (result)
        static_cast<MshModel*>(model())->addMesh(std::move(result));
    else
        OGSError::box("Error adding layer to mesh.");
}

void MshView::extractSurfaceMesh()
{
    QModelIndex const index = this->selectionModel()->currentIndex();
    if (!index.isValid())
        return;

    MeshLib::Mesh const*const mesh = static_cast<MshModel*>(this->model())->getMesh(index);
    SurfaceExtractionDialog dlg;
    if (dlg.exec() != QDialog::Accepted)
        return;

    MathLib::Vector3 const& dir (dlg.getNormal());
    int const tolerance (dlg.getTolerance());
    std::unique_ptr<MeshLib::Mesh> sfc_mesh(
        MeshLib::MeshSurfaceExtraction::getMeshSurface(*mesh, dir, tolerance));
    if (sfc_mesh)
        static_cast<MshModel*>(model())->addMesh(std::move(sfc_mesh));
    else
        OGSError::box(" No surfaces found to extract\n using the specified parameters.");
}

void MshView::convertMeshToGeometry()
{
    QModelIndex const index = this->selectionModel()->currentIndex();
    MeshLib::Mesh const*const mesh = static_cast<MshModel*>(this->model())->getMesh(index);
    emit requestMeshToGeometryConversion(mesh);
}

void MshView::exportToShapefile() const
{
    QModelIndex const index = this->selectionModel()->currentIndex();
    if (!index.isValid())
        return;

    QSettings const settings;
    QFileInfo const fi (settings.value("lastOpenedMeshFileDirectory").toString());
    MeshLib::Mesh const*const mesh = static_cast<MshModel*>(this->model())->getMesh(index);
    QString const fileName = QFileDialog::getSaveFileName(nullptr, "Convert mesh to shapefile...",
        LastSavedFileDirectory::getDir() + QString::fromStdString(mesh->getName()),
        "ESRI Shapefile (*.shp)");
    if (!fileName.isEmpty())
    {
        LastSavedFileDirectory::setDir(fileName);
        if (!FileIO::SHPInterface::write2dMeshToSHP(fileName.toStdString(), *mesh))
            OGSError::box("Error exporting mesh\n to shapefile");
    }
}

void MshView::exportToTetGen()
{
    QModelIndex const index = this->selectionModel()->currentIndex();
    if (!index.isValid())
        return;

    MeshLib::Mesh const*const mesh = static_cast<MshModel*>(this->model())->getMesh(index);
    QSettings const settings;
    QString const filename = QFileDialog::getSaveFileName(this, "Write TetGen input file to",
        settings.value("lastOpenedTetgenFileDirectory").toString(),
        "TetGen Geometry (*.smesh)");
    if (!filename.isEmpty())
    {
        FileIO::TetGenInterface tg;
        std::vector<MeshLib::Node> attr;
        tg.writeTetGenSmesh(filename.toStdString(), *mesh, attr);
    }
}

void MshView::writeToFile() const
{
    QModelIndex const index = this->selectionModel()->currentIndex();
    if (!index.isValid())
    {
        OGSError::box("No mesh selected.");
        return;
    }

    MeshLib::Mesh const*const mesh = static_cast<MshModel*>(this->model())->getMesh(index);
    if (mesh == nullptr)
    {
        OGSError::box("No mesh selected.");
        return;
    }

    SaveMeshDialog dlg(*mesh);
    dlg.exec();
}

void MshView::addDIRECTSourceTerms()
{
    //QModelIndex const index = this->selectionModel()->currentIndex();
    //MeshLib::Mesh const*const grid = static_cast<MshModel*>(this->model())->getMesh(index);
}

void MshView::loadDIRECTSourceTerms()
{
    QModelIndex const index = this->selectionModel()->currentIndex();
    emit loadFEMCondFileRequested(index.data(0).toString().toStdString());
}

void MshView::checkMeshQuality ()
{
    QModelIndex const index = this->selectionModel()->currentIndex();
    MshItem const*const item = static_cast<MshItem*>(static_cast<MshModel*>(this->model())->getItem(index));
    emit qualityCheckRequested(item->vtkSource());
}

