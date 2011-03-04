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
#include <QAbstractItemModel>

//image to mesh conversion
#include "GridAdapter.h"
#include <vtkDataObject.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include "VtkGeoImageSource.h"

#include "msh_mesh.h"
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkGenericDataObjectReader.h>

VtkVisPipelineView::VtkVisPipelineView( QWidget* parent /*= 0*/ )
: QTreeView(parent)
{
	this->setItemsExpandable(false);
	//setEditTriggers(QAbstractItemView::AllEditTriggers);
	CheckboxDelegate* checkboxDelegate = new CheckboxDelegate(this);
	this->setItemDelegateForColumn(1, checkboxDelegate);
	this->header()->setStretchLastSection(false);
	this->header()->setResizeMode(QHeaderView::ResizeToContents);
}

void VtkVisPipelineView::setModel(QAbstractItemModel* model)
{
	QTreeView::setModel(model);
	
	// Move Visisble checkbox to the left.
	// This is done here because at constructor time there arent any sections.
	this->header()->moveSection(1, 0);
}

void VtkVisPipelineView::contextMenuEvent( QContextMenuEvent* event )
{
	QModelIndex index = selectionModel()->currentIndex();
	if (index.isValid())
	{
		// check if object is an image data object
		int objectType = static_cast<VtkVisPipelineItem*>(static_cast<VtkVisPipeline*>(this->model())->getItem(this->selectionModel()->currentIndex()))->algorithm()->GetOutputDataObject(0)->GetDataObjectType();
		bool isSourceItem = (this->selectionModel()->currentIndex().parent().isValid()) ? 0 : 1;

		QMenu menu;
		QAction* addFilterAction = menu.addAction("Add filter...");
		QAction* addMeshingAction = NULL;
		if (objectType == VTK_IMAGE_DATA) 
		{
			isSourceItem = false; // this exception is needed as image object are only displayed in the vis-pipeline
			addMeshingAction = menu.addAction("Convert Image to Mesh...");
			connect(addMeshingAction, SIGNAL(triggered()), this, SLOT(convertImageToMesh()));
		}
		QAction* addConvertToCFEMeshAction = NULL;
		if (objectType == VTK_UNSTRUCTURED_GRID) 
		{
			addConvertToCFEMeshAction = menu.addAction("Convert to Mesh...");
			connect(addConvertToCFEMeshAction, SIGNAL(triggered()), this, SLOT(convertVTKToOGSMesh()));
		}
		menu.addSeparator();
		QAction* exportVtkAction = menu.addAction("Export as VTK");
		QAction* exportOsgAction = menu.addAction("Export as OpenSG");
		QAction* removeAction = NULL;
		if (!isSourceItem) 
		{
			removeAction = menu.addAction("Remove");
			connect(removeAction, SIGNAL(triggered()), this, SLOT(removeSelectedPipelineItem()));
		}

		connect(addFilterAction, SIGNAL(triggered()), this, SLOT(addPipelineFilterItem()));
		connect(exportVtkAction, SIGNAL(triggered()), this, SLOT(exportSelectedPipelineItemAsVtk()));
		connect(exportOsgAction, SIGNAL(triggered()), this, SLOT(exportSelectedPipelineItemAsOsg()));

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

	vtkSmartPointer<VtkGeoImageSource> imageSource = VtkGeoImageSource::SafeDownCast(algorithm);
	vtkSmartPointer<vtkImageData> image = imageSource->GetOutput();

	Mesh_Group::CFEMesh* mesh = GridAdapter::convertImgToMesh(image, imageSource->getOrigin(), imageSource->getSpacing());
	// now do something with the mesh (save, display, whatever... )
	std::string msh_name("NewMesh");
	emit meshAdded(mesh, msh_name);
}

void VtkVisPipelineView::convertVTKToOGSMesh()
{
	vtkSmartPointer<vtkAlgorithm> algorithm = static_cast<VtkVisPipelineItem*>(static_cast<VtkVisPipeline*>(this->model())->getItem(this->selectionModel()->currentIndex()))->algorithm();
	vtkUnstructuredGrid* grid(NULL);

	vtkGenericDataObjectReader* dataReader = vtkGenericDataObjectReader::SafeDownCast(algorithm); // for old filetypes
	if (dataReader) grid = vtkUnstructuredGrid::SafeDownCast(dataReader->GetOutput());
	else
	{
		vtkXMLUnstructuredGridReader* xmlReader = vtkXMLUnstructuredGridReader::SafeDownCast(algorithm); // for new filetypes
		grid = vtkUnstructuredGrid::SafeDownCast(xmlReader->GetOutput());
	}

	Mesh_Group::CFEMesh* mesh = GridAdapter::convertUnstructuredGrid(grid);
	std::string msh_name("NewMesh");
	emit meshAdded(mesh, msh_name);
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
		emit actorSelected(NULL);
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
