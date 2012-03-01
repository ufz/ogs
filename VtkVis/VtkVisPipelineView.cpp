/**
 * \file VtkVisPipelineView.cpp
 * 18/2/2010 LB Initial implementation
 *
 * Implementation of VtkVisPipelineView
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
#include "msh_mesh.h"
#include "GridAdapter.h"
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
		vtkAlgorithm* algorithm =
		        static_cast<VtkVisPipelineItem*>(static_cast<VtkVisPipeline*>(this->model())
		                                         ->
		                                         getItem(this->selectionModel()->
		                                                 currentIndex()))->algorithm();
		int objectType = algorithm->GetOutputDataObject(0)->GetDataObjectType();
		VtkAlgorithmProperties* vtkProps = dynamic_cast<VtkAlgorithmProperties*>(algorithm);
		bool isSourceItem =
		        (this->selectionModel()->currentIndex().parent().isValid()) ? 0 : 1;

		QMenu menu;
		QAction* addFilterAction = menu.addAction("Add filter...");

		QAction* addLUTAction(NULL);
		QAction* addMeshingAction(NULL);
		if (objectType == VTK_IMAGE_DATA)
		{
			isSourceItem = false; // this exception is needed as image object are only displayed in the vis-pipeline
			addMeshingAction = menu.addAction("Convert Image to Mesh...");
			connect(addMeshingAction, SIGNAL(triggered()), this,
			        SLOT(showImageToMeshConversionDialog()));
		}
		else
		{
			addLUTAction = menu.addAction("Add color table...");
			connect(addLUTAction, SIGNAL(triggered()), this, SLOT(addColorTable()));
		}

		QAction* addConvertToCFEMeshAction(NULL);
		if (objectType == VTK_UNSTRUCTURED_GRID)
		{
			addConvertToCFEMeshAction = menu.addAction("Convert to Mesh...");
			connect(addConvertToCFEMeshAction, SIGNAL(triggered()), this,
			        SLOT(convertVTKToOGSMesh()));
		}
		menu.addSeparator();
		QAction* exportVtkAction = menu.addAction("Export as VTK");
		QAction* exportOsgAction = menu.addAction("Export as OpenSG");
		QAction* removeAction = NULL;
		if (!isSourceItem || vtkProps == NULL)
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

		menu.exec(event->globalPos());
	}
}

void VtkVisPipelineView::exportSelectedPipelineItemAsVtk()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
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
	QSettings settings("UFZ", "OpenGeoSys-5");
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
	
	GridAdapter* mesh = VtkMeshConverter::convertImgToMesh(imageSource->GetOutput(), origin, spacing[0], element_type, intensity_type);
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
	GridAdapter* mesh = VtkMeshConverter::convertUnstructuredGrid(grid);
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

	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName = QFileDialog::getOpenFileName(this, "Select color table",
	                                                settings.value("lastOpenedTextureFileDirectory"). toString(), 
													"Color table files (*.xml);;");
	QFileInfo fi(fileName);

	if (fi.suffix().toLower() == "xml")
	{
		VtkVisPointSetItem* pointSetItem = dynamic_cast<VtkVisPointSetItem*>(item);
		if (pointSetItem)
		{
			VtkAlgorithmProperties* props = pointSetItem->getVtkProperties();
			if (props)
			{
				props->SetLookUpTable(array_name, fileName);
				item->SetActiveAttribute(array_name);
				emit requestViewUpdate();
			}
		}
		else
			QMessageBox::warning(NULL, "Color lookup table could not be applied.",
								 "Color lookup tables can only be applied to VtkVisPointSetItem.");
	}
}
