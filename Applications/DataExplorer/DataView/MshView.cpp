/**
 * \file
 * \author Lars Bilke
 * \date   2009-09-24
 * \brief  Implementation of the MshView class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MshView.h"

#include <QHeaderView>
#include <QContextMenuEvent>
#include <QFileDialog>
#include <QMenu>
#include <QObject>
#include <QSettings>

#include "Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLayerEditDialog.h"
#include "MeshValueEditDialog.h"
#include "SurfaceExtractionDialog.h"
#include "MshItem.h"
#include "MshModel.h"
#include "OGSError.h"
#include "MeshSurfaceExtraction.h"

#include "ImportFileTypes.h"
#include "LastSavedFileDirectory.h"
#include "SaveMeshDialog.h"

#include "SHPInterface.h"
#include "TetGenInterface.h"

MshView::MshView( QWidget* parent /*= 0*/ )
	: QTreeView(parent)
{
	setUniformRowHeights(true);
	//resizeColumnsToContents();
	//resizeRowsToContents();
}

MshView::~MshView()
{
}

void MshView::updateView()
{
	setAlternatingRowColors(true);
	setColumnWidth(0,125);
	size_t nColumns = (this->model() != NULL) ? this->model()->columnCount() : 0;
	for (size_t i = 1; i < nColumns; i++)
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
			emit meshSelected(list_item->getMesh());
		}
		else
		{
			emit enableSaveButton(false);
			emit enableRemoveButton(false);
			emit elementSelected(static_cast<InSituLib::VtkMappedMeshSource*>(
				static_cast<const MshItem*>(tree_item->parentItem())->vtkSource()),
				static_cast<unsigned>(tree_item->row()), true);
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
	QModelIndex index = this->selectionModel()->currentIndex();
	MshItem* item = dynamic_cast<MshItem*>(static_cast<TreeItem*>(index.internalPointer()));

	if (item)
	{
		unsigned mesh_dim (item->getMesh()->getDimension());

		QMenu menu;
		QMenu direct_cond_menu("DIRECT Conditions");
		QAction*    editMeshAction = menu.addAction("Edit mesh...");
		QAction*  editValuesAction = menu.addAction("Edit material groups...");
		QAction* meshQualityAction = menu.addAction("Calculate element quality...");
		QAction* surfaceMeshAction (nullptr);
		QAction* tetgenExportAction (nullptr);
		if (mesh_dim==3)
		{
			surfaceMeshAction = menu.addAction("Extract surface...");
			tetgenExportAction = menu.addAction("Export to TetGen...");
		}
		QAction* mesh2geoAction (nullptr);
		QAction* shapeExportAction (nullptr);
		if (mesh_dim==2)
		{
			mesh2geoAction = menu.addAction("Convert to geometry");
			shapeExportAction = menu.addAction("Export to Shapefile...");
		}

		menu.addSeparator();
		//menu.addMenu(&direct_cond_menu);
		QAction*   addDirectAction = direct_cond_menu.addAction("Add...");
		QAction*  loadDirectAction = direct_cond_menu.addAction("Load...");
		//menu.addSeparator();
		connect(editMeshAction,         SIGNAL(triggered()), this, SLOT(openMeshEditDialog()));
		connect(editValuesAction,       SIGNAL(triggered()), this, SLOT(openValuesEditDialog()));
		connect(meshQualityAction,      SIGNAL(triggered()), this, SLOT(checkMeshQuality()));
		if (mesh_dim==3)
		{
			connect(surfaceMeshAction,  SIGNAL(triggered()), this, SLOT(extractSurfaceMesh()));
			connect(tetgenExportAction, SIGNAL(triggered()), this, SLOT(exportToTetGen()));
		}
		connect(addDirectAction,	    SIGNAL(triggered()), this, SLOT(addDIRECTSourceTerms()));
		connect(loadDirectAction,       SIGNAL(triggered()), this, SLOT(loadDIRECTSourceTerms()));
		if (mesh_dim==2)
		{
			connect(mesh2geoAction,     SIGNAL(triggered()), this, SLOT(convertMeshToGeometry()));
			connect(shapeExportAction,  SIGNAL(triggered()), this, SLOT(exportToShapefile()));
		}
		menu.exec(event->globalPos());
	}
}

void MshView::openMeshEditDialog()
{
	MshModel* model = static_cast<MshModel*>(this->model());
	QModelIndex index = this->selectionModel()->currentIndex();
	const MeshLib::Mesh* mesh =
	        static_cast<MshModel*>(this->model())->getMesh(index);

	MeshLayerEditDialog meshLayerEdit(mesh);
	connect(&meshLayerEdit, SIGNAL(mshEditFinished(MeshLib::Mesh*)),
		    model, SLOT(addMesh(MeshLib::Mesh*)));
	meshLayerEdit.exec();
}

void MshView::openValuesEditDialog()
{
	MshModel* model = static_cast<MshModel*>(this->model());
	QModelIndex index = this->selectionModel()->currentIndex();
	MeshLib::Mesh* mesh = const_cast<MeshLib::Mesh*>(static_cast<MshModel*>(this->model())->getMesh(index));

	MeshValueEditDialog valueEdit(mesh);
	connect(&valueEdit, SIGNAL(valueEditFinished(MeshLib::Mesh*)),
		    model, SLOT(updateMesh(MeshLib::Mesh*)));
	valueEdit.exec();
}

void MshView::extractSurfaceMesh()
{
	QModelIndex index = this->selectionModel()->currentIndex();
	if (!index.isValid())
		return;

	const MeshLib::Mesh* mesh = static_cast<MshModel*>(this->model())->getMesh(index);
	SurfaceExtractionDialog dlg;
	if (dlg.exec() != QDialog::Accepted)
		return;
	
	MathLib::Vector3 const& dir (dlg.getNormal());
	int const tolerance (dlg.getTolerance());
	MeshLib::Mesh* sfc_mesh (MeshLib::MeshSurfaceExtraction::getMeshSurface(*mesh, dir, tolerance));
	if (sfc_mesh != nullptr)
		static_cast<MshModel*>(this->model())->addMesh(sfc_mesh);
	else
		OGSError::box(" No surfaces found to extract\n using the specified parameters.");
}

void MshView::convertMeshToGeometry()
{
	QModelIndex index = this->selectionModel()->currentIndex();
	const MeshLib::Mesh* mesh = static_cast<MshModel*>(this->model())->getMesh(index);
	emit requestMeshToGeometryConversion(mesh);
}

void MshView::exportToShapefile() const
{
	QModelIndex index = this->selectionModel()->currentIndex();
	if (!index.isValid())
		return;

	QSettings settings;
	QFileInfo fi (settings.value("lastOpenedMeshFileDirectory").toString());
	const MeshLib::Mesh* mesh = static_cast<MshModel*>(this->model())->getMesh(index);
	QString fileName = QFileDialog::getSaveFileName(NULL, "Convert mesh to shapefile...",
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
	QModelIndex index = this->selectionModel()->currentIndex();

	if (!index.isValid())
		return;

	const MeshLib::Mesh* mesh = static_cast<MshModel*>(this->model())->getMesh(index);
	QSettings settings;
	QString filename = QFileDialog::getSaveFileName(this, "Write TetGen input file to",
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
	QModelIndex index = this->selectionModel()->currentIndex();

	if (!index.isValid())
	{
		OGSError::box("No mesh selected.");
		return;
	}

	const MeshLib::Mesh* mesh = static_cast<MshModel*>(this->model())->getMesh(index);

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
	QModelIndex index = this->selectionModel()->currentIndex();
	const MeshLib::Mesh* grid = static_cast<MshModel*>(this->model())->getMesh(index);
	emit requestCondSetupDialog(grid->getName(), GeoLib::GEOTYPE::INVALID, 0, false);
}

void MshView::loadDIRECTSourceTerms()
{
	QModelIndex index = this->selectionModel()->currentIndex();
	emit loadFEMCondFileRequested(index.data(0).toString().toStdString());
}

void MshView::checkMeshQuality ()
{
	QModelIndex index = this->selectionModel()->currentIndex();
	MshItem* item = static_cast<MshItem*>(static_cast<MshModel*>(this->model())->getItem(index));
	emit qualityCheckRequested(static_cast<InSituLib::VtkMappedMeshSource*>(item->vtkSource()));
}

