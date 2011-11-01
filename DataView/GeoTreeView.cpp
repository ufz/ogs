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

#include "XMLInterface.h"

GeoTreeView::GeoTreeView(QWidget* parent) : QTreeView(parent)
{
//	setContextMenuPolicy(Qt::CustomContextMenu);
//    connect(this, SIGNAL(customContextMenuRequested(const QPoint &)), this, SLOT(showContextMenu(const QPoint &)));
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
	emit itemSelectionChanged(selected, deselected);
	return QTreeView::selectionChanged(selected, deselected);
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
			QAction* addNameAction = menu.addAction("Set name for element...");
			connect(addNameAction, SIGNAL(triggered()), this, SLOT(setNameForElement()));
		}
		// The current index refers to the name of a geometry-object
		else if (item->child(0)->data(0).toString().compare("Points") == 0) // clumsy way to find out
		{
			QAction* addCNDAction = menu.addAction("Add FEM Conditions...");
			QAction* saveAction = menu.addAction("Save geometry...");
			menu.addSeparator();
			QAction* removeAction = menu.addAction("Remove geometry");
			connect(addCNDAction, SIGNAL(triggered()), this, SLOT(addFEMConditions()));
			connect(saveAction, SIGNAL(triggered()), this, SLOT(writeToFile()));
			connect(removeAction, SIGNAL(triggered()), this, SLOT(removeList()));
		}
	}

	menu.exec(event->globalPos());
}

void GeoTreeView::addFEMConditions()
{
	TreeItem* item = static_cast<GeoTreeModel*>(model())->getItem(
	        this->selectionModel()->currentIndex());
	emit loadFEMCondFileRequested(item->data(0).toString().toStdString());
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

void GeoTreeView::setNameForElement()
{
	TreeItem* item = static_cast<GeoTreeModel*>(model())->getItem(
	        this->selectionModel()->currentIndex());
	size_t id = item->data(0).toInt();
	std::string type = GEOLIB::convertGeoTypeToString(dynamic_cast<GeoObjectListItem*>(item->parentItem())->getType());
		std::string geometry_name = item->parentItem()->parentItem()->data(0).toString().toStdString();
	emit requestNameChangeDialog(geometry_name, type, id);
}

void GeoTreeView::writeToFile() const
{
	TreeItem* item = static_cast<GeoTreeModel*>(model())->getItem(
	        this->selectionModel()->currentIndex());
	QString gliName = item->data(0).toString();
	QString fileName = QFileDialog::getSaveFileName(NULL,
						"Save geometry as", gliName, "GeoSys mesh file (*.gml)");
	if (!fileName.isEmpty())
		emit saveToFileRequested(gliName, fileName);
}
