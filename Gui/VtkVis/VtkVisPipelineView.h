/**
 * \file VtkVisPipelineView.h
 * 18/2/2010 LB Initial implementation
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

namespace MeshLib
{
class CFEMesh;
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

	/// Sends an requestRemovePipelineItem() signal to remove
	/// the currently selected item.
	void removeSelectedPipelineItem();

	/// Sends a requestAddPipelineFilterItem() signal to add a filter.
	void addPipelineFilterItem();

	/// Calls the conversion method for creating an OGS Mesh from a vtkImageData object.
	void constructMeshFromImage(QString msh_name, MshElemType::type element_type, UseIntensityAs::type intensity_type);

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
