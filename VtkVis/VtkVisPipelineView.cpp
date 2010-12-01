/**
 * \file VtkVisPipelineView.cpp
 * 18/2/2010 LB Initial implementation
 *
 * Implementation of VtkVisPipelineView
 */

// ** INCLUDES **
#include "VtkVisPipelineView.h"

#include "VtkVisPipelineItem.h"
#include "VtkVisPipeline.h"
#include "CheckboxDelegate.h"

#include <vtkProp3D.h>
#include <vtkDataSetMapper.h>

#include <QMenu>
#include <QContextMenuEvent>
#include <QFileDialog>
#include <QSettings>
#include <QHeaderView>

//image to mesh conversion
#include "GridAdapter.h"
#include <vtkDataObject.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include "VtkGeoImageSource.h"


VtkVisPipelineView::VtkVisPipelineView( QWidget* parent /*= 0*/ )
: QTreeView(parent)
{
	setItemsExpandable(false);
	//setEditTriggers(QAbstractItemView::AllEditTriggers);
	CheckboxDelegate* checkboxDelegate = new CheckboxDelegate(this);
	setItemDelegateForColumn(1, checkboxDelegate);
	header()->setStretchLastSection(false);
	header()->setResizeMode(QHeaderView::ResizeToContents);
}

void VtkVisPipelineView::contextMenuEvent( QContextMenuEvent* event )
{
	QModelIndex index = selectionModel()->currentIndex();
	if (index.isValid())
	{
		// check if object is an image data object
		//int objectType = static_cast<VtkVisPipelineItem*>(static_cast<VtkVisPipeline*>(this->model())->getItem(this->selectionModel()->currentIndex()))->algorithm()->GetOutputDataObject(0)->GetDataObjectType();

		QMenu menu;
		QAction* addFilterAction = menu.addAction("Add filter...");
		//QAction* addMeshingAction = NULL;
		//if (objectType == VTK_IMAGE_DATA) addMeshingAction = menu.addAction("Convert Image to Mesh...");
		menu.addSeparator();
		QAction* exportVtkAction = menu.addAction("Export as VTK");
		QAction* exportOsgAction = menu.addAction("Export as OpenSG");
		QAction* removeAction = menu.addAction("Remove");

		connect(addFilterAction, SIGNAL(triggered()), this, SLOT(addPipelineFilterItem()));
		//connect(addMeshingAction, SIGNAL(triggered()), this, SLOT(convertImageToMesh()));
		connect(exportVtkAction, SIGNAL(triggered()), this, SLOT(exportSelectedPipelineItemAsVtk()));
		connect(exportOsgAction, SIGNAL(triggered()), this, SLOT(exportSelectedPipelineItemAsOsg()));
		connect(removeAction, SIGNAL(triggered()), this, SLOT(removeSelectedPipelineItem()));

		menu.exec(event->globalPos());
	}
}

void VtkVisPipelineView::exportSelectedPipelineItemAsVtk()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QModelIndex idx = this->selectionModel()->currentIndex();
	QString filename = QFileDialog::getSaveFileName(this, "Export object to vtk-file",
		settings.value("lastExportedFileDirectory").toString(),"VTK file (*.vtk)");
	if (!filename.isEmpty())
	{
		static_cast<VtkVisPipelineItem*>(static_cast<VtkVisPipeline*>(this->model())->getItem(idx))->writeToFile(filename.toStdString());
		QDir dir = QDir(filename);
		settings.setValue("lastExportedFileDirectory", dir.absolutePath());
	}
}

void VtkVisPipelineView::exportSelectedPipelineItemAsOsg()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QModelIndex idx = this->selectionModel()->currentIndex();
	QString filename = QFileDialog::getSaveFileName(this, "Export object to OpenSG file",
		settings.value("lastExportedFileDirectory").toString(), "OpenSG file (*.osb)");
	if (!filename.isEmpty())
	{
		static_cast<VtkVisPipelineItem*>(static_cast<VtkVisPipeline*>(this->model())->getItem(idx))->writeToFile(filename.toStdString());
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

void VtkVisPipelineView::convertImageToMesh()
{
	vtkSmartPointer<vtkAlgorithm> algorithm = static_cast<VtkVisPipelineItem*>(static_cast<VtkVisPipeline*>(this->model())->getItem(this->selectionModel()->currentIndex()))->algorithm();

	vtkSmartPointer<VtkGeoImageSource> imageSource = vtkSmartPointer<VtkGeoImageSource>(VtkGeoImageSource::SafeDownCast(algorithm));
	vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>(imageSource->GetOutput());

	Mesh_Group::CFEMesh* mesh = GridAdapter::convertImgToMesh(image, imageSource->getOrigin(), imageSource->getSpacing());
	// now do something with the mesh (save, display, whatever... )
}

void VtkVisPipelineView::selectionChanged( const QItemSelection &selected, const QItemSelection &deselected )
{
	QTreeView::selectionChanged(selected, deselected);

	QModelIndex index = this->selectionModel()->currentIndex();
	if (index.isValid())
	{
		VtkVisPipelineItem* item = static_cast<VtkVisPipelineItem*>(index.internalPointer());
		emit actorSelected(item->actor());
		emit itemSelected(item);
	}
	else
	{
		emit actorSelected((vtkProp3D*)NULL);
		emit itemSelected(NULL);
	}
		
}

void VtkVisPipelineView::selectItem( vtkProp3D* actor )
{
	QModelIndex index = ((VtkVisPipeline*)(this->model()))->getIndex(actor);
	if (!index.isValid())
		return;

	blockSignals(true);
	QItemSelectionModel* selectionModel = this->selectionModel();
	selectionModel->clearSelection();
	selectionModel->select(index, QItemSelectionModel::Select);
	blockSignals(false);
}
