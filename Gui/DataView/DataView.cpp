/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file DataView.cpp
 *
 * Created on 2009-09-24 by Lars Bilke
 */

#include "DataView.h"
#include "Mesh.h"
#include "MshEditDialog.h"
#include "MshItem.h"
#include "MshModel.h"
#include "OGSError.h"
#include <QHeaderView>

#include "VtkMeshSource.h"

#include <QContextMenuEvent>
#include <QFileDialog>
#include <QMenu>
#include <QObject>
#include <QSettings>

#include "Legacy/MeshIO.h"
#include "XmlIO/VTKInterface.h"
#include "Writer.h" // necessary to avoid Linker Error in Windows

DataView::DataView( QWidget* parent /*= 0*/ )
	: QTreeView(parent)
{
	//resizeColumnsToContents();
	//resizeRowsToContents();
}

DataView::~DataView()
{
}

void DataView::updateView()
{
	setAlternatingRowColors(true);
	setColumnWidth(0,125);
	size_t nColumns = (this->model() != NULL) ? this->model()->columnCount() : 0;
	for (size_t i = 1; i < nColumns; i++)
		resizeColumnToContents(i);
}

void DataView::addMeshAction()
{
	QSettings settings;
	QString fileName =
	        QFileDialog::getOpenFileName(this, "Select mesh file", settings.value(
	                                             "lastOpenedFileDirectory").toString(),
	                                     "OpenGeosys mesh files (*.msh);;All files (* *.*)");
	if (!fileName.isEmpty())
	{
		std::string name = fileName.toStdString();
		FileIO::MeshIO meshIO;
		MeshLib::Mesh* msh = meshIO.loadMeshFromFile(name);
		if (msh)
			static_cast<MshModel*>(this->model())->addMesh(msh);
	}
}

void DataView::removeMesh()
{
	emit requestMeshRemoval(this->selectionModel()->currentIndex());
}

void DataView::removeAllMeshes()
{
	TreeItem* root = static_cast<MshModel*>(this->model())->getItem(QModelIndex());
	int nChildren = root->childCount() - 1;
	for (int i = nChildren; i >= 0; i--)
		emit requestMeshRemoval(this->model()->index(i, 0, QModelIndex()));
}

void DataView::contextMenuEvent( QContextMenuEvent* event )
{
	QModelIndex index = this->selectionModel()->currentIndex();
	MshItem* item = dynamic_cast<MshItem*>(static_cast<TreeItem*>(index.internalPointer()));

	if (item)
	{
		QMenu menu;
		QAction* editMeshAction   = menu.addAction("Edit mesh...");
		QAction* checkMeshAction  = menu.addAction("Check mesh quality...");
		QAction* saveMeshAction   = menu.addAction("Save mesh...");
		menu.addSeparator();
		QMenu direct_cond_menu("DIRECT Conditions");
		menu.addMenu(&direct_cond_menu);
		QAction* addDirectAction  = direct_cond_menu.addAction("Add...");
		QAction* loadDirectAction = direct_cond_menu.addAction("Load...");
		menu.addSeparator();
		QAction* removeMeshAction = menu.addAction("Remove mesh");
		connect(editMeshAction, SIGNAL(triggered()), this, SLOT(openMshEditDialog()));
		connect(checkMeshAction, SIGNAL(triggered()), this, SLOT(checkMeshQuality()));
		connect(saveMeshAction, SIGNAL(triggered()), this, SLOT(writeMeshToFile()));
		connect(addDirectAction, SIGNAL(triggered()), this, SLOT(addDIRECTSourceTerms()));
		connect(loadDirectAction, SIGNAL(triggered()), this, SLOT(loadDIRECTSourceTerms()));
		connect(removeMeshAction, SIGNAL(triggered()), this, SLOT(removeMesh()));
		menu.exec(event->globalPos());
	}
}

void DataView::openMshEditDialog()
{
	MshModel* model = static_cast<MshModel*>(this->model());
	QModelIndex index = this->selectionModel()->currentIndex();
	const MeshLib::Mesh* mesh =
	        static_cast<MshModel*>(this->model())->getMesh(index);

	MshEditDialog meshEdit(mesh);
	connect(&meshEdit, SIGNAL(mshEditFinished(MeshLib::Mesh*)),
		    model, SLOT(addMesh(MeshLib::Mesh*)));
	meshEdit.exec();
}

int DataView::writeMeshToFile() const
{
	QModelIndex index = this->selectionModel()->currentIndex();
	const MeshLib::Mesh* mesh =
	        static_cast<MshModel*>(this->model())->getMesh(index);

	if (mesh)
	{
		QSettings settings;
		QString mshName = QString::fromStdString(
		        static_cast<MshModel*>(this->model())->getMesh(index)->getName());
		QString fileName = QFileDialog::getSaveFileName(NULL, "Save mesh as",
		                                    settings.value("lastOpenedMeshFileDirectory").toString(),
											"VTK Unstructured Grid (*.vtu);;GeoSys legacy mesh file (*.msh)");

		if (!fileName.isEmpty())
		{
			QFileInfo fi(fileName);
			if (fi.suffix().toLower() == "vtu")
			{
				FileIO::VTKInterface vtkIO;
				vtkIO.setMesh(mesh);
				vtkIO.writeToFile(fileName.toStdString().c_str());
			}
			if (fi.suffix().toLower() == "msh")
			{
				FileIO::MeshIO meshIO;
				meshIO.setMesh(mesh);
				meshIO.writeToFile(fileName.toStdString().c_str());
			}
			QDir dir = QDir(fileName);
			settings.setValue("lastOpenedMeshFileDirectory", dir.absolutePath());
			return 1;
		}
		else
			OGSError::box("No file name entered.");
	}
	return 0;
}

void DataView::addDIRECTSourceTerms()
{
	QModelIndex index = this->selectionModel()->currentIndex();
	const MeshLib::Mesh* grid = static_cast<MshModel*>(this->model())->getMesh(index);
	emit requestCondSetupDialog(grid->getName(), GeoLib::INVALID, 0, false);
}


void DataView::loadDIRECTSourceTerms()
{
	QModelIndex index = this->selectionModel()->currentIndex();
	const MeshLib::Mesh* grid = static_cast<MshModel*>(this->model())->getMesh(index);
	// TODO6 const std::vector<MeshLib::Node*>* nodes = grid->getNodes();
	// TODO6 emit requestDIRECTSourceTerms(grid->getName(), nodes);
}

void DataView::checkMeshQuality ()
{
	QModelIndex index = this->selectionModel()->currentIndex();
	MshItem* item = static_cast<MshItem*>(static_cast<MshModel*>(this->model())->getItem(index));
	emit qualityCheckRequested(item->vtkSource());
}

/*
   void DataView::selectionChanged( const QItemSelection &selected, const QItemSelection &deselected )
   {
    emit itemSelectionChanged(selected, deselected);
    return QTreeView::selectionChanged(selected, deselected);
   }

   void DataView::selectionChangedFromOutside( const QItemSelection &selected, const QItemSelection &deselected )
   {
    QItemSelectionModel* selModel = this->selectionModel();

    Q_ASSERT(selModel);

    selModel->blockSignals(true);
    selModel->select(deselected, QItemSelectionModel::Deselect);
    selModel->select(selected, QItemSelectionModel::Select);
    selModel->blockSignals(false);

    Model* model = static_cast<Model*>(this->model());
    //model->setSelectionFromOutside(selected, deselected);

    return QTreeView::selectionChanged(selected, deselected);
   }

   void DataView::clearSelection()
   {
    selectionModel()->clearSelection();
   }
 */

