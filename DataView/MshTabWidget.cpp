/**
 * \file MshTabWidget.cpp
 * 3/11/2009 LB Initial implementation
 * 18/05/2010 KR Re-Implementation
 *
 * Implementation of MshTabWidget
 */

// ** INCLUDES **
#include "GridAdapter.h"
#include "MshEditDialog.h"
#include "MshTabWidget.h"
#include "MshModel.h"
#include "OGSError.h"
#include "TreeItem.h"
#include "MeshQualityChecker.h"

#include <QObject>
#include <QMenu>
#include <QContextMenuEvent>
#include <QFileDialog>
#include <QSettings>
#include <QInputDialog>

MshTabWidget::MshTabWidget( QWidget* parent /*= 0*/ )
: QWidget(parent)
{
	setupUi(this);

	connect(this->addMeshPushButton, SIGNAL(clicked()), this, SLOT(addMeshAction()));
	connect(this->clearAllPushButton, SIGNAL(clicked()), this, SLOT(removeAllMeshes()));


/*
	mshTableView->setSelectionMode(QAbstractItemView::ExtendedSelection);
	mshTableView->setSelectionBehavior(QAbstractItemView::SelectRows);
	mshNodeTableView->setSelectionMode(QAbstractItemView::ExtendedSelection);
	mshNodeTableView->setSelectionBehavior(QAbstractItemView::SelectRows);

	connect(mshTableView, SIGNAL(itemSelectionChanged(QItemSelection, QItemSelection)),
		this, SLOT(changeMshSubmodelViews(QItemSelection, QItemSelection)));
*/
}

void MshTabWidget::addMeshAction()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName = QFileDialog::getOpenFileName(this, "Select mesh file", settings.value("lastOpenedFileDirectory").toString(), "OpenGeosys mesh files (*.msh);;All files (* *.*)");
	if (!fileName.isEmpty())
	{
		std::string name = fileName.toStdString();
		Mesh_Group::CFEMesh* msh = MshModel::loadMeshFromFile(name);
		if (msh) static_cast<MshModel*>(this->treeView->model())->addMesh(msh, name);
	}
}

void MshTabWidget::removeMesh()
{
	emit requestMeshRemoval(this->treeView->selectionModel()->currentIndex());
}

void MshTabWidget::removeAllMeshes()
{
	TreeItem* root = static_cast<MshModel*>(this->treeView->model())->getItem(QModelIndex());
	int nChildren = root->childCount()-1;
	for (int i=nChildren; i>=0; i--)
		emit requestMeshRemoval(this->treeView->model()->index(i, 0, QModelIndex()));
}


void MshTabWidget::contextMenuEvent( QContextMenuEvent* event )
{
	QMenu menu;
	QAction* editMeshAction    = menu.addAction("Edit mesh...");
	QAction* checkMeshAction    = menu.addAction("Check mesh quality...");
	QAction* saveMeshAction    = menu.addAction("Save mesh...");
	menu.addSeparator();
	QAction* removeMeshAction    = menu.addAction("Remove...");
	connect(editMeshAction, SIGNAL(triggered()), this, SLOT(openMshEditDialog()));
	connect(checkMeshAction, SIGNAL(triggered()), this, SLOT(checkMeshQuality()));
	connect(saveMeshAction, SIGNAL(triggered()), this, SLOT(writeMeshToFile()));
	connect(removeMeshAction, SIGNAL(triggered()), this, SLOT(removeMesh()));
	menu.exec(event->globalPos());
}

void MshTabWidget::openMshEditDialog()
{
	MshModel* model = static_cast<MshModel*>(this->treeView->model());
	QModelIndex index = this->treeView->selectionModel()->currentIndex();
	const Mesh_Group::CFEMesh* mesh = static_cast<MshModel*>(this->treeView->model())->getMesh(index)->getCFEMesh();

	MshEditDialog meshEdit(mesh);
	connect(&meshEdit, SIGNAL(mshEditFinished(Mesh_Group::CFEMesh*, std::string&)), model, SLOT(addMesh(Mesh_Group::CFEMesh*, std::string&)));
	meshEdit.exec();
}

int MshTabWidget::writeMeshToFile() const
{
	QModelIndex index = this->treeView->selectionModel()->currentIndex();
	const Mesh_Group::CFEMesh* mesh = static_cast<MshModel*>(this->treeView->model())->getMesh(index)->getCFEMesh();

	if (mesh)
	{
		QString mshName = QString::fromStdString(static_cast<MshModel*>(this->treeView->model())->getMesh(index)->getName());
		std::string fileName = QFileDialog::getSaveFileName(NULL, "Save mesh as", mshName, "GeoSys mesh file (*.msh)").toStdString();

		if (!fileName.empty())
		{
			std::fstream* out = new std::fstream(fileName.c_str(), std::fstream::out);
			if (out->is_open())
			{
				mesh->Write(out);
				out->close();
				return 1;
			}
			else
				std::cout << "MshTabWidget::saveMeshFile() - Could not create file..." << std::endl;
		}
		else OGSError::box("No file name entered.");
	}
	return 0;
}

void MshTabWidget::checkMeshQuality () const
{
	QModelIndex index = this->treeView->selectionModel()->currentIndex();
	const Mesh_Group::CFEMesh* mesh = static_cast<MshModel*>(this->treeView->model())->getMesh(index)->getCFEMesh();

	if (mesh) {
		QString msh_name = QString::fromStdString(static_cast<MshModel*>(this->treeView->model())->getMesh(index)->getName());
		std::string file_name = QFileDialog::getSaveFileName(NULL, "Save histogram as", msh_name, "raw data (*.txt)").toStdString();
		Mesh_Group::MeshQualityChecker checker (mesh);
		checker.check ();

		// simple suggestion: number of classes with Sturges criterion
		size_t nclasses (static_cast<size_t>(1 + 3.3 * log ((mesh->getElementVector()).size())));
		bool ok;
		size_t size (static_cast<size_t>(QInputDialog::getInt(NULL, "OGS-Histogramm", "number of histogramm classes/spins (min: 1, max: 10000)", static_cast<int>(nclasses), 1, 10000, 1, &ok)));

		std::vector<size_t> histogramm (size,0);
		checker.getHistogramm(histogramm);
		std::ofstream out (file_name.c_str());
		const size_t histogramm_size (histogramm.size());
		for (size_t k(0); k<histogramm_size; k++) {
			out << k/static_cast<double>(histogramm_size) << " " << histogramm[k] << std::endl;
		}
		out.close ();
	}
}

/*
void MshTabWidget::changeMshSubmodelViews( QItemSelection selected, QItemSelection deselected )
{

	if (selected.size() > 0)
	{
		QModelIndex index = *(selected.begin()->indexes().begin());
		if (!index.isValid())
			return;

		MshModel* mshModel = static_cast<MshModel*>(mshTableView->model());

		ModelItem* item = mshModel->itemFromIndex(index);

		mshNodeTableView->setModel(item->models()[0]);
		mshNodeTableView->resizeColumnsToContents();
		mshNodeTableView->resizeRowsToContents();

		connect(mshNodeTableView, SIGNAL(itemSelectionChanged(const QItemSelection&,const QItemSelection&)),
			item->models()[0], SLOT(setSelection(const QItemSelection&, const QItemSelection&)));

//		connect(item->models()[0], SIGNAL(dataChanged(QModelIndex,QModelIndex)), _scene, SLOT(updateItems(QModelIndex,QModelIndex)));
//		connect(item->models()[1], SIGNAL(dataChanged(QModelIndex,QModelIndex)), _scene, SLOT(updateItems(QModelIndex,QModelIndex)));

		mshElemTableView->setModel(item->models()[1]);
		mshElemTableView->resizeColumnsToContents();
		mshElemTableView->resizeRowsToContents();

	}
	else
	{
		mshNodeTableView->setModel(NULL);
		mshElemTableView->setModel(NULL);
	}

}
*/
