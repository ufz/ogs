/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-18
 * \brief  Implementation of the VtkVisPipelineView class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkVisPipelineView.h"

#include <vtkDataSetMapper.h>
#include <vtkProp3D.h>

#include <QAbstractItemModel>
#include <QContextMenuEvent>
#include <QFileDialog>
#include <QHeaderView>
#include <QMenu>
#include <QMessageBox>
#include <QSettings>

#include "Base/CheckboxDelegate.h"
#include "Base/OGSError.h"
#include "MeshLib/IO/VtkIO/VtkMeshConverter.h"
#include "MeshLib/Mesh.h"
#include "MeshToolsLib/MeshGenerators/RasterToMesh.h"
#include "VtkVisImageItem.h"
#include "VtkVisPipeline.h"
#include "VtkVisPipelineItem.h"
#include "VtkVisPointSetItem.h"

// image to mesh conversion
#include <vtkDataObject.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkTransformFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridAlgorithm.h>
#include <vtkXMLUnstructuredGridReader.h>

#include "MeshFromRasterDialog.h"
#include "VtkGeoImageSource.h"

VtkVisPipelineView::VtkVisPipelineView(QWidget* parent /*= 0*/)
    : QTreeView(parent)
{
    this->setItemsExpandable(false);
    auto* checkboxDelegate = new CheckboxDelegate(this);
    this->setItemDelegateForColumn(1, checkboxDelegate);
    this->header()->setStretchLastSection(true);
    this->header()->setSectionResizeMode(QHeaderView::ResizeToContents);
}

void VtkVisPipelineView::setModel(QAbstractItemModel* model)
{
    QTreeView::setModel(model);

    // Move Visible checkbox to the left.
    // This is done here because at constructor time there aren't any sections.
    this->header()->moveSection(1, 0);
}

void VtkVisPipelineView::contextMenuEvent(QContextMenuEvent* event)
{
    QModelIndex index = selectionModel()->currentIndex();
    if (index.isValid())
    {
        // check object type
        VtkVisPipelineItem* item = static_cast<VtkVisPipelineItem*>(
            static_cast<VtkVisPipeline*>(this->model())
                ->getItem(this->selectionModel()->currentIndex()));
        int objectType =
            item->algorithm()->GetOutputDataObject(0)->GetDataObjectType();
        VtkAlgorithmProperties* vtkProps = item->getVtkProperties();
        bool isSourceItem =
            !(this->selectionModel()->currentIndex().parent().isValid());

        QMenu menu;
        QAction* addFilterAction = menu.addAction("Add filter...");

        if (objectType == VTK_IMAGE_DATA)
        {
            // this exception is needed as image object are only displayed in
            // the vis-pipeline
            isSourceItem = false;
            QAction* saveRasterAction = menu.addAction("Save Raster...");
            connect(saveRasterAction, SIGNAL(triggered()), this,
                    SLOT(writeRaster()));
            QAction* addMeshingAction =
                menu.addAction("Convert Raster to Mesh...");
            connect(addMeshingAction, SIGNAL(triggered()), this,
                    SLOT(showImageToMeshConversionDialog()));
        }
        else
        {
            QAction* addLUTAction = menu.addAction("Add color table...");
            connect(addLUTAction, SIGNAL(triggered()), this,
                    SLOT(addColorTable()));
        }

        if (objectType == VTK_UNSTRUCTURED_GRID)
        {
            QAction* addConvertToMeshAction =
                menu.addAction("Convert to Mesh...");
            connect(addConvertToMeshAction, SIGNAL(triggered()), this,
                    SLOT(convertVTKToOGSMesh()));
        }
        menu.addSeparator();
        QAction* exportVtkAction = menu.addAction("Export as VTK");

        if (!isSourceItem || vtkProps->IsRemovable())
        {
            menu.addSeparator();
            QAction* removeAction = menu.addAction("Remove");
            connect(removeAction, SIGNAL(triggered()), this,
                    SLOT(removeSelectedPipelineItem()));
        }

        connect(addFilterAction, SIGNAL(triggered()), this,
                SLOT(addPipelineFilterItem()));
        connect(exportVtkAction, SIGNAL(triggered()), this,
                SLOT(exportSelectedPipelineItemAsVtk()));

        menu.exec(event->globalPos());
    }
}

void VtkVisPipelineView::exportSelectedPipelineItemAsVtk()
{
    QSettings settings;
    QModelIndex idx = this->selectionModel()->currentIndex();
    QString filename = QFileDialog::getSaveFileName(
        this, "Export object to vtk-file",
        settings.value("lastExportedFileDirectory").toString(),
        "VTK file (*.*)");
    if (!filename.isEmpty())
    {
        static_cast<VtkVisPipelineItem*>(
            static_cast<VtkVisPipeline*>(this->model())->getItem(idx))
            ->writeToFile(filename.toStdString());
        QDir dir = QDir(filename);
        settings.setValue("lastExportedFileDirectory", dir.absolutePath());
    }
}

void VtkVisPipelineView::removeSelectedPipelineItem()
{
    emit requestRemovePipelineItem(selectionModel()->currentIndex());
}

void VtkVisPipelineView::addPipelineFilterItem()
{
    emit requestAddPipelineFilterItem(selectionModel()->currentIndex());
}

void VtkVisPipelineView::showImageToMeshConversionDialog()
{
    MeshFromRasterDialog dlg;
    if (dlg.exec() != QDialog::Accepted)
    {
        return;
    }

    vtkSmartPointer<vtkAlgorithm> algorithm =
        static_cast<VtkVisPipelineItem*>(
            static_cast<VtkVisPipeline*>(this->model())
                ->getItem(this->selectionModel()->currentIndex()))
            ->algorithm();

    vtkSmartPointer<VtkGeoImageSource> imageSource =
        VtkGeoImageSource::SafeDownCast(algorithm);

    auto const raster = VtkGeoImageSource::convertToRaster(imageSource);
    if (!raster)
    {
        OGSError::box("Image could not be converted to a raster.");
        return;
    }

    auto mesh = MeshToolsLib::RasterToMesh::convert(
        *raster, dlg.getElementSelection(), dlg.getIntensitySelection(),
        dlg.getArrayName());
    if (mesh)
    {
        mesh->setName(dlg.getMeshName());
        emit meshAdded(mesh.release());
    }
    else
    {
        OGSError::box("Error creating mesh.");
    }
}

void VtkVisPipelineView::writeRaster()
{
    QModelIndex const index = this->selectionModel()->currentIndex();
    static_cast<VtkVisImageItem*>(
        static_cast<VtkVisPipeline*>(this->model())->getItem(index))
        ->writeAsRaster();
}

void VtkVisPipelineView::convertVTKToOGSMesh()
{
    VtkVisPipelineItem* item = static_cast<VtkVisPipelineItem*>(
        static_cast<VtkVisPipeline*>(this->model())
            ->getItem(this->selectionModel()->currentIndex()));
    vtkSmartPointer<vtkAlgorithm> algorithm = item->algorithm();

    vtkUnstructuredGrid* grid(nullptr);
    vtkUnstructuredGridAlgorithm* ugAlg =
        vtkUnstructuredGridAlgorithm::SafeDownCast(algorithm);
    if (ugAlg)
    {
        grid = ugAlg->GetOutput();
    }
    else
    {
        // for old filetypes
        vtkGenericDataObjectReader* dataReader =
            vtkGenericDataObjectReader::SafeDownCast(algorithm);
        if (dataReader)
        {
            grid = vtkUnstructuredGrid::SafeDownCast(dataReader->GetOutput());
        }
        else
        {
            // for new filetypes
            vtkXMLUnstructuredGridReader* xmlReader =
                vtkXMLUnstructuredGridReader::SafeDownCast(algorithm);
            grid = vtkUnstructuredGrid::SafeDownCast(xmlReader->GetOutput());
        }
    }
    MeshLib::Mesh* mesh =
        MeshLib::VtkMeshConverter::convertUnstructuredGrid(grid);
    mesh->setName(item->data(0).toString().toStdString());
    emit meshAdded(mesh);
}

void VtkVisPipelineView::selectionChanged(const QItemSelection& selected,
                                          const QItemSelection& deselected)
{
    QTreeView::selectionChanged(selected, deselected);

    if (selected.empty())
    {
        return;
    }

    QModelIndex index = *selected.indexes().begin();
    if (index.isValid())
    {
        auto* item = static_cast<VtkVisPipelineItem*>(index.internalPointer());
        emit actorSelected(item->actor());
        emit itemSelected(item);
        if (item->transformFilter())
        {
            emit dataObjectSelected(vtkDataObject::SafeDownCast(
                item->transformFilter()->GetOutputDataObject(0)));
        }
    }
    else
    {
        emit actorSelected(nullptr);
        emit itemSelected(nullptr);
        emit dataObjectSelected(nullptr);
    }
}

void VtkVisPipelineView::selectItem(vtkProp3D* actor)
{
    this->selectItem(
        static_cast<VtkVisPipeline*>(this->model())->getIndex(actor));
}

void VtkVisPipelineView::selectItem(const QModelIndex& index)
{
    if (!index.isValid())
    {
        return;
    }

    QItemSelectionModel* selectionModel = this->selectionModel();
    selectionModel->clearSelection();
    selectionModel->select(index, QItemSelectionModel::Select);
}

void VtkVisPipelineView::addColorTable()
{
    VtkVisPipelineItem* item(static_cast<VtkVisPipelineItem*>(
        static_cast<VtkVisPipeline*>(this->model())
            ->getItem(this->selectionModel()->currentIndex())));
    const QString array_name = item->GetActiveAttribute();

    QSettings settings;
    QString filename = QFileDialog::getOpenFileName(
        this, "Select color table",
        settings.value("lastOpenedLutFileDirectory").toString(),
        "Color table files (*.xml);;");
    QFileInfo fi(filename);

    if (fi.suffix().toLower() == "xml")
    {
        auto* pointSetItem = dynamic_cast<VtkVisPointSetItem*>(item);
        if (pointSetItem)
        {
            VtkAlgorithmProperties* props = pointSetItem->getVtkProperties();
            if (props)
            {
                props->SetLookUpTable(array_name, filename);
                item->SetActiveAttribute(array_name);
                emit requestViewUpdate();
            }
        }
        else
        {
            QMessageBox::warning(nullptr,
                                 "Color lookup table could not be applied.",
                                 "Color lookup tables can only be applied to "
                                 "VtkVisPointSetItem.");
        }
        QDir dir = QDir(filename);
        settings.setValue("lastOpenedLutFileDirectory", dir.absolutePath());
    }
}
