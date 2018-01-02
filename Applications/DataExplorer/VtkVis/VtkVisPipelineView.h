/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-18
 * \brief  Definition of the VtkVisPipelineView class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

// ** INCLUDES **
#include <QTreeView>

class QItemSelection;
class QAbstractItemModel;
class VtkVisPipelineItem;
class vtkProp3D;
class vtkDataObject;

enum class MeshElemType;

namespace MeshLib {
    class Mesh;
    enum class UseIntensityAs;
    enum class MeshElemType;
}

/**
 * \brief VtkVisPipelineView is a QTreeView and shows VtkVisPipelineItems and their relation to each other.
 */
class VtkVisPipelineView : public QTreeView
{
    Q_OBJECT

public:
    /// @brief Constructor.
    VtkVisPipelineView(QWidget* parent = nullptr);

    /// @brief Overridden to set model specific header properties.
    void setModel(QAbstractItemModel* model) override;

protected slots:
    /// Emits itemSelected() signals when an items was selected.
    void selectionChanged(const QItemSelection& selected,
                          const QItemSelection& deselected) override;
    void selectItem(vtkProp3D* actor);
    void selectItem(const QModelIndex &index);

private:
    /// Creates a menu on right-clicking on an item.
    void contextMenuEvent(QContextMenuEvent* event) override;

private slots:
    /// Adds a color lookup table to the current scalar array of the selected pipeline item.
    void addColorTable();

    /// Exports the currently selected item as a VTK file
    void exportSelectedPipelineItemAsVtk();

    /// Exports the currently selected item as a Fbx file.
    void exportSelectedPipelineItemAsFbx();

    /// Sends an requestRemovePipelineItem() signal to remove
    /// the currently selected item.
    void removeSelectedPipelineItem();

    /// Sends a requestAddPipelineFilterItem() signal to add a filter.
    void addPipelineFilterItem();

    /// Calls the dialog to
    void showImageToMeshConversionDialog();

    /// Calls the conversion method for making a vtk grid an ogs mesh.
    void convertVTKToOGSMesh();

signals:
    void requestViewUpdate();
    void requestRemovePipelineItem(QModelIndex);
    void requestAddPipelineFilterItem(QModelIndex);
    void itemSelected(VtkVisPipelineItem*);
    void actorSelected(vtkProp3D*);
    void dataObjectSelected(vtkDataObject*);
    void meshAdded(MeshLib::Mesh*);
};
