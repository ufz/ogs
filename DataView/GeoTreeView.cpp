/**
 * \file GeoTreeView.cpp
 * 2011/02/07 KR Initial implementation
 */

#include <QFileDialog>
#include <QMenu>
#include <iostream>

#include "GeoObjectListItem.h"
#include "GeoTreeItem.h"
#include "GeoTreeModel.h"
#include "GeoTreeView.h"
#include "OGSError.h"

GeoTreeView::GeoTreeView(QWidget* parent) : QTreeView(parent)
{
}

void GeoTreeView::updateView()
{
	setAlternatingRowColors(true);
	//resizeColumnToContents(0);
	setColumnWidth(1,150);
	setColumnWidth(1,50);
	setColumnWidth(2,50);
}

void GeoTreeView::on_Clicked(QModelIndex idx)
{
	qDebug("%d, %d",idx.parent().row(), idx.row());
}

void GeoTreeView::selectionChanged( const QItemSelection &selected,
                                    const QItemSelection &deselected )
{
	Q_UNUSED(deselected);
	if (!selected.isEmpty())
	{
		const QModelIndex idx = *(selected.indexes().begin());
		const TreeItem* tree_item = static_cast<TreeModel*>(this->model())->getItem(idx);

		const GeoObjectListItem* list_item = dynamic_cast<GeoObjectListItem*>(tree_item->parentItem());
		if (list_item)
			emit geoItemSelected(list_item->vtkSource(), tree_item->row());
		else
			emit removeGeoItemSelection();
	}
	//emit itemSelectionChanged(selected, deselected);
	//return QTreeView::selectionChanged(selected, deselected);
}

void GeoTreeView::selectionChangedFromOutside( const QItemSelection &selected,
                                               const QItemSelection &deselected )
{
	QItemSelectionModel* selModel = this->selectionModel();

	selModel->blockSignals(true);
	selModel->select(deselected, QItemSelectionModel::Deselect);
	selModel->select(selected, QItemSelectionModel::Select);
	selModel->blockSignals(false);

	return QTreeView::selectionChanged(selected, deselected);
}

void GeoTreeView::contextMenuEvent( QContextMenuEvent* event )
{
	QModelIndex index = this->selectionModel()->currentIndex();
	TreeItem* item = static_cast<TreeItem*>(index.internalPointer());

	GeoObjectListItem* list = dynamic_cast<GeoObjectListItem*>(item);
	QMenu menu;

	// The current index is a list of points/polylines/surfaces
	if (list != NULL)
	{
		QAction* connectPlyAction(NULL);
		if (list->getType() == GEOLIB::POLYLINE)
		{
			connectPlyAction = menu.addAction("Connect Polylines...");
			connect(connectPlyAction, SIGNAL(triggered()), this,
			        SLOT(connectPolylines()));
		}
		menu.addSeparator();
		QAction* removeAction = menu.addAction("Remove " + item->data(0).toString());
		connect(removeAction, SIGNAL(triggered()), this, SLOT(removeList()));
	}
	else
	{
		GeoObjectListItem* parent = dynamic_cast<GeoObjectListItem*>(item->parentItem());

		// The current index refers to a geo-object
		if (parent != NULL)
		{
			QAction* addCondAction = menu.addAction("Set as FEM condition...");
			QAction* addNameAction = menu.addAction("Set name...");
			connect(addCondAction, SIGNAL(triggered()), this, SLOT(setElementAsCondition()));
			connect(addNameAction, SIGNAL(triggered()), this, SLOT(setNameForElement()));
		}
		// The current index refers to the name of a geometry-object
		else if (item->childCount() > 0)
		{
			if (item->child(0)->data(0).toString().compare("Points") == 0) // clumsy way to find out
			{
				QAction* saveAction = menu.addAction("Save geometry...");
				QAction* addCNDAction = menu.addAction("Load FEM Conditions...");
				//QAction* saveCondAction    = menu.addAction("Save FEM conditions...");
				menu.addSeparator();
				QAction* removeAction = menu.addAction("Remove geometry");
				connect(saveAction, SIGNAL(triggered()), this, SLOT(writeToFile()));
				connect(addCNDAction, SIGNAL(triggered()), this, SLOT(loadFEMConditions()));
				//connect(saveCondAction, SIGNAL(triggered()), this, SLOT(saveFEMConditions()));
				connect(removeAction, SIGNAL(triggered()), this, SLOT(removeList()));
			}
		}
	}

	menu.exec(event->globalPos());
}

void GeoTreeView::connectPolylines()
{
	TreeItem* item = static_cast<GeoTreeModel*>(model())->getItem(
	        this->selectionModel()->currentIndex())->parentItem();
	emit requestLineEditDialog(item->data(0).toString().toStdString());
}

void GeoTreeView::removeList()
{
	TreeItem* item = static_cast<GeoTreeModel*>(model())->getItem(
	        this->selectionModel()->currentIndex());

	GeoObjectListItem* list = dynamic_cast<GeoObjectListItem*>(item);
	if (list)
		emit listRemoved((item->parentItem()->data(
		                          0).toString()).toStdString(), list->getType());
	else
		emit listRemoved((item->data(0).toString()).toStdString(), GEOLIB::INVALID);
}

void GeoTreeView::setElementAsCondition()
{
	const TreeItem* item = static_cast<GeoTreeModel*>(model())->getItem(
	        this->selectionModel()->currentIndex());
	const size_t id = item->row();
	const GEOLIB::GEOTYPE type = static_cast<GeoObjectListItem*>(item->parentItem())->getType();
	const std::string geometry_name = item->parentItem()->parentItem()->data(0).toString().toStdString();
	emit requestCondSetupDialog(geometry_name, type, id);
}

void GeoTreeView::setNameForElement()
{
	const TreeItem* item = static_cast<GeoTreeModel*>(model())->getItem(
	        this->selectionModel()->currentIndex());
	const size_t id = item->row();
	const GEOLIB::GEOTYPE type = static_cast<GeoObjectListItem*>(item->parentItem())->getType();
	const std::string geometry_name = item->parentItem()->parentItem()->data(0).toString().toStdString();
	emit requestNameChangeDialog(geometry_name, type, id);
}

void GeoTreeView::writeToFile() const
{
	TreeItem* item = static_cast<GeoTreeModel*>(model())->getItem(
	        this->selectionModel()->currentIndex());
	QString gliName = item->data(0).toString();
	QString fileName = QFileDialog::getSaveFileName(NULL,
						"Save geometry as", gliName, "GeoSys geometry file (*.gml)");
	if (!fileName.isEmpty())
		emit saveToFileRequested(gliName, fileName);
}

void GeoTreeView::loadFEMConditions()
{
	TreeItem* item = static_cast<GeoTreeModel*>(model())->getItem(
	        this->selectionModel()->currentIndex());
	emit loadFEMCondFileRequested(item->data(0).toString().toStdString());
}
/*
void GeoTreeView::saveFEMConditions()
{
	TreeItem* item = static_cast<GeoTreeModel*>(model())->getItem(
	        this->selectionModel()->currentIndex());
	QString fileName = QFileDialog::getSaveFileName(NULL,
						"Save FEM Conditions as", "", "OpenGeoSys FEM Condition file (*.cnd);; GeoSys Boundary Condition (*.bc);; GeoSys Initial Condition (*.ic);; GeoSys Source Condition (*.st)");
	emit saveFEMConditionsRequested(item->data(0).toString(), fileName);
}
*/