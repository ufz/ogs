/**
 * \file GeoTreeView.cpp
 * 2011/02/07 KR Initial implementation
 */

#include <iostream>
#include <QFileDialog>
#include <QMenu>

#include "GeoTreeView.h"
#include "GeoTreeModel.h"
#include "GeoTreeItem.h"
#include "GeoObjectListItem.h"
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

void GeoTreeView::selectionChanged( const QItemSelection &selected, const QItemSelection &deselected )
{
	emit itemSelectionChanged(selected, deselected);
	return QTreeView::selectionChanged(selected, deselected);
}

void GeoTreeView::selectionChangedFromOutside( const QItemSelection &selected, const QItemSelection &deselected )
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

	if (item->childCount() > 0)
	{
		// The current index refers to the name of a geometry-object
		if (item->child(0)->data(0).toString().compare("Points") == 0) // clumsy way to find out
		{
			QAction* addCNDAction = menu.addAction("Add FEM Conditions...");
			QAction* saveAction = menu.addAction("Save geometry...");
			menu.addSeparator();
			QAction* removeAction = menu.addAction("Remove geometry");
			connect(addCNDAction, SIGNAL(triggered()), this, SLOT(addFEMConditions()));
			connect(saveAction, SIGNAL(triggered()), this, SLOT(writeToFile()));
			connect(removeAction, SIGNAL(triggered()), this, SLOT(removeList()));
		}
		else if (list != NULL)
		{
			QAction* connectPlyAction(NULL);
			if (list->getType() == GEOLIB::POLYLINE)
			{
				connectPlyAction = menu.addAction("Connect Polylines...");
				connect(connectPlyAction, SIGNAL(triggered()), this, SLOT(connectPolylines()));
			}
			menu.addSeparator();
			QAction* removeAction = menu.addAction("Remove " + item->data(0).toString());
			connect(removeAction, SIGNAL(triggered()), this, SLOT(removeList()));
		}
		// The current index refers to a geo object
		else
		{
			QString temp_name;
			QMenu menu;
	/*
			if (static_cast<GeoTreeModel*>(model())->objectFromIndex(index, temp_name)->type() == GEOLIB::POINT)
			{
				QAction* stratAction = menu.addAction("Display Stratigraphy...");
				QAction* exportAction = menu.addAction("Export to GMS...");
				connect(stratAction, SIGNAL(triggered()), this, SLOT(displayStratigraphy()));
				connect(exportAction, SIGNAL(triggered()), this, SLOT(exportStation()));
				menu.exec(event->globalPos());
			}
			else
			{
				menu.addAction("View Information...");
				QAction* showDiagramAction = menu.addAction("View Diagram...");
				connect(showDiagramAction, SIGNAL(triggered()), this, SLOT(showDiagramPrefsDialog()));
				menu.exec(event->globalPos());
			}
	*/
		}
	}
	menu.exec(event->globalPos());
}

void GeoTreeView::addFEMConditions()
{
	TreeItem* item = static_cast<GeoTreeModel*>(model())->getItem(this->selectionModel()->currentIndex());
	emit loadFEMCondFileRequested(item->data(0).toString().toStdString());
}

void GeoTreeView::writeToFile() const
{
	TreeItem* item = static_cast<GeoTreeModel*>(model())->getItem(this->selectionModel()->currentIndex());
	QString gliName = item->data(0).toString();
	QString fileName = QFileDialog::getSaveFileName(NULL, "Save geometry as", gliName, "GeoSys mesh file (*.gml)");
	if (!fileName.isEmpty()) {
		emit saveToFileRequested(gliName, fileName);
	}
}

void GeoTreeView::removeList()
{
	TreeItem* item = static_cast<GeoTreeModel*>(model())->getItem(this->selectionModel()->currentIndex());

	GeoObjectListItem* list = dynamic_cast<GeoObjectListItem*>(item);
	if (list) emit listRemoved((item->parentItem()->data(0).toString()).toStdString(), list->getType());
	else emit listRemoved((item->data(0).toString()).toStdString(), GEOLIB::INVALID);
}

void GeoTreeView::connectPolylines()
{
	TreeItem* item = static_cast<GeoTreeModel*>(model())->getItem(this->selectionModel()->currentIndex())->parentItem();
	emit requestLineEditDialog(item->data(0).toString().toStdString());
}

