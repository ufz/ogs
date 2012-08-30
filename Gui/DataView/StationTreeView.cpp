/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file StationTreeView.cpp
 *
 * Created on by Karsten Rink
 */

#include <QFileDialog>
#include <QMenu>
#include <iostream>

//TODO6 #include "GMSInterface.h"
#include "Station.h"
//TODO6 #include "StationIO.h"

#include "DiagramPrefsDialog.h"
#include "ModelTreeItem.h"
#include "OGSError.h"
#include "StationTreeModel.h"
#include "StationTreeView.h"
#include "StratWindow.h"

StationTreeView::StationTreeView(QWidget* parent) : QTreeView(parent)
{
//	setContextMenuPolicy(Qt::CustomContextMenu);
//    connect(this, SIGNAL(customContextMenuRequested(const QPoint &)), this, SLOT(showContextMenu(const QPoint &)));
}

void StationTreeView::updateView()
{
	setAlternatingRowColors(true);
	resizeColumnToContents(0);
	setColumnWidth(1,50);
	setColumnWidth(2,50);
}

void StationTreeView::on_Clicked(QModelIndex idx)
{
	qDebug("%d, %d",idx.parent().row(), idx.row());
}

void StationTreeView::selectionChanged( const QItemSelection &selected,
                                        const QItemSelection &deselected )
{
	emit itemSelectionChanged(selected, deselected);
	return QTreeView::selectionChanged(selected, deselected);
}

void StationTreeView::selectionChangedFromOutside( const QItemSelection &selected,
                                                   const QItemSelection &deselected )
{
	QItemSelectionModel* selModel = this->selectionModel();

	selModel->blockSignals(true);
	selModel->select(deselected, QItemSelectionModel::Deselect);
	selModel->select(selected, QItemSelectionModel::Select);
	selModel->blockSignals(false);

	return QTreeView::selectionChanged(selected, deselected);
}

void StationTreeView::contextMenuEvent( QContextMenuEvent* event )
{
	QModelIndex index = this->selectionModel()->currentIndex();
	ModelTreeItem* item = static_cast<ModelTreeItem*>(index.internalPointer());

	if (!item)  // Otherwise sometimes it crashes when (unmotivated ;-) ) clicking in a treeview
		return;

	// The current index refers to a parent item (e.g. a listname)
	if (item->childCount() > 0)
	{
		QMenu menu;
		QAction* propertyAction = menu.addAction("Display list properties...");
		QAction* exportAction   = menu.addAction("Export to GMS...");
		QAction* saveAction   = menu.addAction("Save to file...");
		menu.addSeparator();
		QAction* removeAction   = menu.addAction("Remove station list");

		connect(propertyAction, SIGNAL(triggered()), this, SLOT(showPropertiesDialog()));
		connect(exportAction,   SIGNAL(triggered()), this, SLOT(exportList()));
		connect(saveAction,   SIGNAL(triggered()), this, SLOT(saveList()));
		connect(removeAction,   SIGNAL(triggered()), this, SLOT(removeStationList()));
		menu.exec(event->globalPos());
	}
	// The current index refers to a station object
	else
	{
		QString temp_name;
		QMenu menu;

		if (static_cast<StationTreeModel*>(model())->stationFromIndex(index,
		                                                              temp_name)->type() ==
		    GeoLib::Station::BOREHOLE)
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
			connect(showDiagramAction, SIGNAL(triggered()), this,
			        SLOT(showDiagramPrefsDialog()));
			menu.exec(event->globalPos());
		}
	}
}

void StationTreeView::displayStratigraphy()
{
	QModelIndex index = this->selectionModel()->currentIndex();

	QString temp_name;
	// get list name
	static_cast<StationTreeModel*>(model())->stationFromIndex(
	        this->selectionModel()->currentIndex(), temp_name);
	// get color table (horrible way to do it but there you go ...)
	std::map<std::string,
	         GeoLib::Color*> colorLookupTable =
	        static_cast<VtkStationSource*>(static_cast<StationTreeModel*>(model())->vtkSource(
	                                               temp_name.
	                                               toStdString()))
	        ->getColorLookupTable();
	StratWindow* stratView =
	        new StratWindow(static_cast<GeoLib::StationBorehole*>(static_cast<StationTreeModel*>(
	                                                                      model())
	                                                              ->stationFromIndex(index,
	                                                                                 temp_name)),
	                        &colorLookupTable);
	stratView->setAttribute(Qt::WA_DeleteOnClose); // this fixes the memory leak shown by cppcheck
	stratView->show();
}

void StationTreeView::saveList()
{
	TreeItem* item = static_cast<StationTreeModel*>(model())->getItem(
	        this->selectionModel()->currentIndex());
	QString listName = item->data(0).toString();
	QString fileName = QFileDialog::getSaveFileName(this, "Save station list", "","*.stn");
	if (!fileName.isEmpty())
		emit stationListSaved(listName, fileName);
}

void StationTreeView::exportList()
{
	// only a test for the stratigraphy screenshot tool and writer!!
	//QString Name = static_cast<StationTreeModel*>(model())->getItem(this->selectionModel()->currentIndex())->data(0).toString();
	//writeStratigraphiesAsImages(Name);

	TreeItem* item = static_cast<StationTreeModel*>(model())->getItem(
	        this->selectionModel()->currentIndex());
	std::string listName = item->data(0).toString().toStdString();
	QString fileName = QFileDialog::getSaveFileName(this,
	                                                "Export Boreholes to GMS-Format",
	                                                "",
	                                                "*.txt");
	if (!fileName.isEmpty())
		emit stationListExportRequested(listName, fileName.toStdString());
}

void StationTreeView::exportStation()
{
	QModelIndex index = this->selectionModel()->currentIndex();
	QString fileName = QFileDialog::getSaveFileName(this,
	                                                "Export Borehole to GMS-Format",
	                                                "",
	                                                "*.txt");
	if (!fileName.isEmpty())
	{
		QString temp_name;
		std::vector<std::string> temp_soil_names;
		temp_soil_names.push_back(""); // soil name vector needs to be initialised
		/* TODO6
		GMSInterface::writeBoreholeToGMS(static_cast<GeoLib::StationBorehole*>(static_cast<
		                                                                               StationTreeModel
		                                                                               *>(
		                                                                               model())->stationFromIndex(index,
		                                                                                                          temp_name)),
		                                 fileName.toStdString(), temp_soil_names);
		*/
	}
}

void StationTreeView::removeStationList()
{
	TreeItem* item = static_cast<StationTreeModel*>(model())->getItem(
	        this->selectionModel()->currentIndex());
	emit stationListRemoved((item->data(0).toString()).toStdString());
}

void StationTreeView::showPropertiesDialog()
{
	QModelIndex index = this->selectionModel()->currentIndex();
	QString name = (static_cast<ModelTreeItem*>(index.internalPointer())->data(0)).toString();
	emit propertiesDialogRequested(name.toStdString());
}

void StationTreeView::showDiagramPrefsDialog()
{
	QModelIndex index = this->selectionModel()->currentIndex();
	emit diagramRequested(index);
}

void StationTreeView::writeStratigraphiesAsImages(QString listName)
{
	std::map<std::string,
	         GeoLib::Color*> colorLookupTable =
	        static_cast<VtkStationSource*>(static_cast<StationTreeModel*>(model())->vtkSource(
	                                               listName.
	                                               toStdString()))
	        ->getColorLookupTable();
	std::vector<ModelTreeItem*> lists = static_cast<StationTreeModel*>(model())->getLists();
	size_t nLists = lists.size();
	for (size_t i = 0; i < nLists; i++)
		if ( listName.toStdString().compare( lists[i]->data(0).toString().toStdString() )
		     == 0 )
		{
			const std::vector<GeoLib::Point*>* stations =
			        dynamic_cast<BaseItem*>(lists[i]->getItem())->getStations();

			for (size_t i = 0; i < stations->size(); i++)
			{
				StratWindow* stratView =
				        new StratWindow(static_cast<GeoLib::StationBorehole*>((*
				                                                               stations)
				                                                              [i]),
				                        &colorLookupTable);
				stratView->setAttribute(Qt::WA_DeleteOnClose); // this fixes the memory leak shown by cppcheck
				stratView->show();
				stratView->stationView->saveAsImage(
				        "c:/project/" +
				        QString::fromStdString(static_cast<GeoLib::StationBorehole*>((
				                                                                             *
				                                                                             stations)
				                                                                     [
				                                                                             i])->getName()) + ".jpg");
				stratView->close();
			}
		}
}
