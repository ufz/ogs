/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-18
 * \brief  Definition of the VtkVisPipelineView class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VTKVISPIPELINEVIEW_H
#define VTKVISPIPELINEVIEW_H

// ** INCLUDES **
#include <QTreeView>
#include "VtkMeshConverter.h"

class QItemSelection;
class QAbstractItemModel;
class VtkVisPipelineItem;
class vtkProp3D;
class vtkDataObject;

namespace MeshLib {
	class Mesh;
}

/**
 * \brief VtkVisPipelineView is a QTreeView and shows VtkVisPipelineItems and their relation to each other.
 */
class VtkVisPipelineView : public QTreeView
{
	Q_OBJECT

public:
	/// @brief Constructor.
	VtkVisPipelineView(QWidget* parent = 0);

	/// @brief Overridden to set model specific header properties.
	virtual void setModel(QAbstractItemModel* model);

protected slots:
	/// Emits itemSelected() signals when an items was selected.
	void selectionChanged(const QItemSelection &selected, const QItemSelection &deselected);
	void selectItem(vtkProp3D* actor);

private:
	/// Creates a menu on right-clicking on an item.
	void contextMenuEvent(QContextMenuEvent* event);

private slots:
	/// Adds a color lookup table to the current scalar array of the selected pipeline item.
	void addColorTable();

	/// Exports the currently selected item as a VTK file
	void exportSelectedPipelineItemAsVtk();

	/// Exports the currently selected item as an OpenSG file
	void exportSelectedPipelineItemAsOsg();

#ifdef VTKFBXCONVERTER_FOUND
	/// Exports the currently selected item as a Fbx file.
	void exportSelectedPipelineItemAsFbx();
#endif // VTKFBXCONVERTER_FOUND

	/// Sends an requestRemovePipelineItem() signal to remove
	/// the currently selected item.
	void removeSelectedPipelineItem();

	/// Sends a requestAddPipelineFilterItem() signal to add a filter.
	void addPipelineFilterItem();

	/// Calls the conversion method for creating an OGS Mesh from a vtkImageData object.
	void constructMeshFromImage(QString msh_name, MeshElemType element_type, UseIntensityAs intensity_type);

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

#endif // VTKVISPIPELINEVIEW_H
