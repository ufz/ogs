/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-18
 * \brief  Implementation of the VtkVisPipelineView class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkVisPipelineView.h"

#include "CheckboxDelegate.h"
#include "VtkVisPipeline.h"
#include "VtkVisPipelineItem.h"
#include "VtkVisPointSetItem.h"

#include <vtkDataSetMapper.h>
#include <vtkProp3D.h>

#include <QAbstractItemModel>
#include <QContextMenuEvent>
#include <QFileDialog>
#include <QHeaderView>
#include <QMenu>
#include <QSettings>
#include <QMessageBox>

//image to mesh conversion
#include "Mesh.h"
#include "VtkGeoImageSource.h"
#include <vtkImageData.h>
#include "MeshFromRasterDialog.h"
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkTransformFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridAlgorithm.h>
#include <vtkXMLUnstructuredGridReader.h>


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
		// check object type
		VtkVisPipelineItem* item = static_cast<VtkVisPipelineItem*>(static_cast<VtkVisPipeline*>(
			this->model())->getItem(this->selectionModel()->currentIndex()));
		int objectType = item->algorithm()->GetOutputDataObject(0)->GetDataObjectType();
		VtkAlgorithmProperties* vtkProps = item->getVtkProperties();
		bool isSourceItem =
		        (this->selectionModel()->currentIndex().parent().isValid()) ? 0 : 1;

		QMenu menu;
		QAction* addFilterAction = menu.addAction("Add filter...");

		QAction* addLUTAction(NULL);
		QAction* addMeshingAction(NULL);
		if (objectType == VTK_IMAGE_DATA)
		{
			// this exception is needed as image object are only displayed in the vis-pipeline
			isSourceItem = false;
			addMeshingAction = menu.addAction("Convert Image to Mesh...");
			connect(addMeshingAction, SIGNAL(triggered()), this,
			        SLOT(showImageToMeshConversionDialog()));
		}
		else
		{
			addLUTAction = menu.addAction("Add color table...");
			connect(addLUTAction, SIGNAL(triggered()), this, SLOT(addColorTable()));
		}

		QAction* addConvertToMeshAction(NULL);
		if (objectType == VTK_UNSTRUCTURED_GRID)
		{
			addConvertToMeshAction = menu.addAction("Convert to Mesh...");
			connect(addConvertToMeshAction, SIGNAL(triggered()), this,
			        SLOT(convertVTKToOGSMesh()));
		}
		menu.addSeparator();
		QAction* exportVtkAction = menu.addAction("Export as VTK");
		QAction* exportOsgAction = menu.addAction("Export as OpenSG");
#ifdef VTKFBXCONVERTER_FOUND
		QAction* exportFbxAction = menu.addAction("Export as Fbx");
#endif
		QAction* removeAction = NULL;
		if (!isSourceItem || vtkProps->IsRemovable())
		{
			removeAction = menu.addAction("Remove");
			connect(removeAction, SIGNAL(triggered()), this,
			        SLOT(removeSelectedPipelineItem()));
		}

		connect(addFilterAction, SIGNAL(triggered()), this, SLOT(addPipelineFilterItem()));
		connect(exportVtkAction, SIGNAL(triggered()), this,
		        SLOT(exportSelectedPipelineItemAsVtk()));
		connect(exportOsgAction, SIGNAL(triggered()), this,
		        SLOT(exportSelectedPipelineItemAsOsg()));
#ifdef VTKFBXCONVERTER_FOUND
		connect(exportFbxAction, SIGNAL(triggered()), this,
		        SLOT(exportSelectedPipelineItemAsFbx()));
#endif

		menu.exec(event->globalPos());
	}
}

void VtkVisPipelineView::exportSelectedPipelineItemAsVtk()
{
	QSettings settings;
	QModelIndex idx = this->selectionModel()->currentIndex();
	QString filename = QFileDialog::getSaveFileName(this, "Export object to vtk-file",
	                                settings.value("lastExportedFileDirectory").toString(),
									"All files (* *.*)");
	if (!filename.isEmpty())
	{
		static_cast<VtkVisPipelineItem*>(static_cast<VtkVisPipeline*>(this->model())->
		                                 getItem(idx))->writeToFile(filename.toStdString());
		QDir dir = QDir(filename);
		settings.setValue("lastExportedFileDirectory", dir.absolutePath());
	}
}

void VtkVisPipelineView::exportSelectedPipelineItemAsOsg()
{
	QSettings settings;
	QModelIndex idx = this->selectionModel()->currentIndex();
	QString filename = QFileDialog::getSaveFileName(this, "Export object to OpenSG file",
	                                                settings.value("lastExportedFileDirectory").
	                                                toString(), "OpenSG file (*.osb)");
	if (!filename.isEmpty())
	{
		static_cast<VtkVisPipelineItem*>(static_cast<VtkVisPipeline*>(this->model())->
		                                 getItem(idx))->writeToFile(filename.toStdString());
		QDir dir = QDir(filename);
		settings.setValue("lastExportedFileDirectory", dir.absolutePath());
	}
}

#ifdef VTKFBXCONVERTER_FOUND
void VtkVisPipelineView::exportSelectedPipelineItemAsFbx()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QModelIndex idx = this->selectionModel()->currentIndex();
	QString filename = QFileDialog::getSaveFileName(this, "Export object to Fbx file",
	                                                settings.value("lastExportedFileDirectory").
	                                                toString(), "Fbx file (*.fbx)");
	if (!filename.isEmpty())
	{
		static_cast<VtkVisPipelineItem*>(static_cast<VtkVisPipeline*>(this->model())->
		                                 getItem(idx))->writeToFile(filename.toStdString());
		QDir dir = QDir(filename);
		settings.setValue("lastExportedFileDirectory", dir.absolutePath());
	}
}
#endif // VTKFBXCONVERTER_FOUND

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
	MeshFromRasterDialog* dlg = new MeshFromRasterDialog();
	connect(dlg, SIGNAL(setMeshParameters(QString, MshElemType::type, UseIntensityAs::type)),
		    this, SLOT(constructMeshFromImage(QString, MshElemType::type, UseIntensityAs::type)));
	dlg->exec();
}

void VtkVisPipelineView::constructMeshFromImage(QString msh_name, MshElemType::type element_type, UseIntensityAs::type intensity_type)
{
	vtkSmartPointer<vtkAlgorithm> algorithm =
	        static_cast<VtkVisPipelineItem*>(static_cast<VtkVisPipeline*>(this->model())->
	                                         getItem(this->selectionModel()->currentIndex()))->algorithm();

	vtkSmartPointer<VtkGeoImageSource> imageSource = VtkGeoImageSource::SafeDownCast(algorithm);
	double origin[3];
	imageSource->GetOutput()->GetOrigin(origin);
	double spacing[3];
	imageSource->GetOutput()->GetSpacing(spacing);

	MeshLib::Mesh* mesh = VtkMeshConverter::convertImgToMesh(imageSource->GetOutput(), origin, spacing[0], element_type, intensity_type);
	mesh->setName(msh_name.toStdString());
	emit meshAdded(mesh);
}

void VtkVisPipelineView::convertVTKToOGSMesh()
{
	VtkVisPipelineItem* item = static_cast<VtkVisPipelineItem*>(static_cast<VtkVisPipeline*>(this->model())->getItem(
												this->selectionModel()->currentIndex()));
	vtkSmartPointer<vtkAlgorithm> algorithm = item->algorithm();


	vtkUnstructuredGrid* grid(NULL);
	vtkUnstructuredGridAlgorithm* ugAlg = vtkUnstructuredGridAlgorithm::SafeDownCast(algorithm);
	if (ugAlg)
		grid = ugAlg->GetOutput();
	else
	{
		// for old filetypes
		vtkGenericDataObjectReader* dataReader = vtkGenericDataObjectReader::SafeDownCast(algorithm);
		if (dataReader)
			grid = vtkUnstructuredGrid::SafeDownCast(dataReader->GetOutput());
		else
		{
			// for new filetypes
			vtkXMLUnstructuredGridReader* xmlReader = vtkXMLUnstructuredGridReader::SafeDownCast(algorithm);
			grid = vtkUnstructuredGrid::SafeDownCast(xmlReader->GetOutput());
		}
	}
	MeshLib::Mesh* mesh = VtkMeshConverter::convertUnstructuredGrid(grid);
	mesh->setName(item->data(0).toString().toStdString());
	emit meshAdded(mesh);
}

void VtkVisPipelineView::selectionChanged( const QItemSelection &selected,
                                           const QItemSelection &deselected )
{
	QTreeView::selectionChanged(selected, deselected);

	QModelIndex index = this->selectionModel()->currentIndex();
	if (index.isValid())
	{
		VtkVisPipelineItem* item = static_cast<VtkVisPipelineItem*>(index.internalPointer());
		emit actorSelected(item->actor());
		emit itemSelected(item);
		if (item->transformFilter())
			emit dataObjectSelected(vtkDataObject::SafeDownCast(
			                                item->transformFilter()->
			                                GetOutputDataObject(0)));
	}
	else
	{
		emit actorSelected(NULL);
		emit itemSelected(NULL);
		emit dataObjectSelected(NULL);
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

void VtkVisPipelineView::addColorTable()
{
	VtkVisPipelineItem* item (
	        static_cast<VtkVisPipelineItem*>(static_cast<VtkVisPipeline*>(this->model())->
	                                         getItem(
	                                                 this->selectionModel()->currentIndex())) );
	const QString array_name = item->GetActiveAttribute();

	QSettings settings;
	QString filename = QFileDialog::getOpenFileName(this, "Select color table",
	                                                settings.value("lastOpenedLutFileDirectory"). toString(),
													"Color table files (*.xml);;");
	QFileInfo fi(filename);

	if (fi.suffix().toLower() == "xml")
	{
		VtkVisPointSetItem* pointSetItem = dynamic_cast<VtkVisPointSetItem*>(item);
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
			QMessageBox::warning(NULL, "Color lookup table could not be applied.",
								 "Color lookup tables can only be applied to VtkVisPointSetItem.");
		QDir dir = QDir(filename);
		settings.setValue("lastOpenedLutFileDirectory", dir.absolutePath());
	}
}
