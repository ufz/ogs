/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Implementation of the StationTreeView class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <QFileDialog>
#include <QMenu>

#include "GMSInterface.h"
#include "Station.h"

#include "DiagramPrefsDialog.h"
#include "ModelTreeItem.h"
#include "OGSError.h"
#include "StationTreeModel.h"
#include "StationTreeView.h"
#include "StratWindow.h"
#include "ImportFileTypes.h"
#include "LastSavedFileDirectory.h"

StationTreeView::StationTreeView(QWidget* parent) : QTreeView(parent)
{
//    setContextMenuPolicy(Qt::CustomContextMenu);
//    connect(this, SIGNAL(customContextMenuRequested(const QPoint &)), this, SLOT(showContextMenu(const QPoint &)));
}

void StationTreeView::updateView()
{
    setAlternatingRowColors(true);
    resizeColumnToContents(0);
    setColumnWidth(1,50);
    setColumnWidth(2,50);
}

void StationTreeView::selectionChanged( const QItemSelection &selected,
                                        const QItemSelection &deselected )
{
    Q_UNUSED(deselected);
    if (!selected.isEmpty())
    {
        const QModelIndex idx = *(selected.indexes().begin());
        const TreeItem* tree_item = static_cast<TreeModel*>(this->model())->getItem(idx);

        const auto* list_item =
            dynamic_cast<const ModelTreeItem*>(tree_item->parentItem());
        if (list_item->getItem())
        {
            if (list_item)
                emit geoItemSelected(list_item->getItem()->vtkSource(), tree_item->row());
            emit enableRemoveButton(false);
            emit enableSaveButton(false);
        }
        else
        {
            emit removeGeoItemSelection();
            emit enableSaveButton(true);
            emit enableRemoveButton(true);
        }
    }
    //emit itemSelectionChanged(selected, deselected);
    //return QTreeView::selectionChanged(selected, deselected);
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
    auto* item = static_cast<ModelTreeItem*>(index.internalPointer());

    if (!item)  // Otherwise sometimes it crashes when (unmotivated ;-) ) clicking in a treeview
        return;

    // The current index refers to a parent item (e.g. a listname)
    if (item->childCount() > 0)
    {
        QMenu menu;
        QAction* exportAction   = menu.addAction("Export to GMS...");
        menu.addSeparator();

        connect(exportAction,   SIGNAL(triggered()), this, SLOT(exportList()));
        menu.exec(event->globalPos());
    }
    // The current index refers to a station object
    else
    {
        QString temp_name;
        QMenu menu;

        QAction* setNameAction = menu.addAction("Set name...");
        connect(setNameAction, SIGNAL(triggered()), this, SLOT(setNameForElement()));
        if (static_cast<StationTreeModel*>(model())->stationFromIndex(index, temp_name)->type()
                == GeoLib::Station::StationType::BOREHOLE)
        {
            QAction* stratAction = menu.addAction("Display Stratigraphy...");
            QAction* exportAction = menu.addAction("Export to GMS...");
            connect(stratAction, SIGNAL(triggered()), this, SLOT(displayStratigraphy()));
            connect(exportAction, SIGNAL(triggered()), this, SLOT(exportStation()));
            menu.exec(event->globalPos());
        }
        else
        {
            QAction* showDiagramAction = menu.addAction("View Diagram...");
            connect(showDiagramAction, SIGNAL(triggered()), this,
                    SLOT(showDiagramPrefsDialog()));
            menu.exec(event->globalPos());
        }
    }
}

void StationTreeView::setNameForElement()
{
    TreeItem const*const item = static_cast<StationTreeModel*>(model())->getItem(
            this->selectionModel()->currentIndex());
    std::string const stn_vec_name = item->parentItem()->data(0).toString().toStdString();
    emit requestNameChangeDialog(stn_vec_name, item->row());
}

void StationTreeView::displayStratigraphy()
{
    QModelIndex index = this->selectionModel()->currentIndex();

    QString temp_name;
    // get list name
    static_cast<StationTreeModel*>(model())->stationFromIndex(
            this->selectionModel()->currentIndex(), temp_name);
    // get color table (horrible way to do it but there you go ...)
    std::map<std::string, DataHolderLib::Color> colorLookupTable =
        static_cast<VtkStationSource*>(static_cast<StationTreeModel*>
            (model())->vtkSource(temp_name.toStdString()))->getColorLookupTable();
    auto* stratView = new StratWindow(
        static_cast<GeoLib::StationBorehole*>(
            static_cast<StationTreeModel*>(model())->stationFromIndex(
                index, temp_name)),
        &colorLookupTable);
    stratView->setAttribute(Qt::WA_DeleteOnClose); // this fixes the memory leak shown by cppcheck
    stratView->show();
}

void StationTreeView::addStationList()
{
    emit openStationListFile(ImportFileType::OGS_STN);
}

void StationTreeView::writeToFile()
{
    QModelIndex index (this->selectionModel()->currentIndex());
    if (!index.isValid())
        OGSError::box("No station list selected.");
    else
    {
        TreeItem* item = static_cast<StationTreeModel*>(model())->getItem(index);
        QString listName = item->data(0).toString();
        QString fileName = QFileDialog::getSaveFileName(this, "Save station list",
            LastSavedFileDirectory::getDir() + listName, "*.stn");
        if (!fileName.isEmpty())
        {
            LastSavedFileDirectory::setDir(fileName);
            emit stationListSaved(listName, fileName);
        }
    }
}

void StationTreeView::exportList()
{
    // only a test for the stratigraphy screenshot tool and writer!!
    //QString Name = static_cast<StationTreeModel*>(model())->getItem(this->selectionModel()->currentIndex())->data(0).toString();
    //writeStratigraphiesAsImages(Name);

    TreeItem* item = static_cast<StationTreeModel*>(model())->getItem(
            this->selectionModel()->currentIndex());
    QString listName = item->data(0).toString();
    QString fileName = QFileDialog::getSaveFileName(this, "Export Boreholes to GMS-Format",
        LastSavedFileDirectory::getDir() + listName, "*.txt");
    if (!fileName.isEmpty())
    {
        LastSavedFileDirectory::setDir(fileName);
        emit stationListExportRequested(listName.toStdString(), fileName.toStdString());
    }
}

void StationTreeView::exportStation()
{
    QModelIndex index = this->selectionModel()->currentIndex();
    QString fileName = QFileDialog::getSaveFileName(this, "Export Borehole to GMS-Format",
        LastSavedFileDirectory::getDir(), "*.txt");
    if (!fileName.isEmpty())
    {
        QString temp_name;
        std::vector<GeoLib::Point*> stations;
        stations.push_back(static_cast<GeoLib::StationBorehole*>(
                    static_cast<StationTreeModel*>(model())->stationFromIndex(index,temp_name)));
        FileIO::GMSInterface::writeBoreholesToGMS(&stations, fileName.toStdString());
        LastSavedFileDirectory::setDir(fileName);
    }
}

void StationTreeView::removeStationList()
{
    QModelIndex index (this->selectionModel()->currentIndex());
    if (!index.isValid())
        OGSError::box("No station list selected.");
    else
    {
        TreeItem* item = static_cast<StationTreeModel*>(model())->getItem(index);
        emit stationListRemoved((item->data(0).toString()).toStdString());

        if(this->selectionModel()->selectedIndexes().count() == 0)
        {
            emit enableSaveButton(false);
            emit enableRemoveButton(false);
        }
    }
}

void StationTreeView::showDiagramPrefsDialog()
{
    QModelIndex index = this->selectionModel()->currentIndex();
    emit diagramRequested(index);
}

void StationTreeView::writeStratigraphiesAsImages(QString listName)
{
    std::map<std::string, DataHolderLib::Color> colorLookupTable =
        static_cast<VtkStationSource*>(static_cast<StationTreeModel*>
            (model())->vtkSource(listName.toStdString()))->getColorLookupTable();
    std::vector<ModelTreeItem*> lists = static_cast<StationTreeModel*>(model())->getLists();
    std::size_t nLists = lists.size();
    for (std::size_t i = 0; i < nLists; i++)
    {
        if ( listName.compare( lists[i]->data(0).toString() ) != 0 )
            continue;

        std::vector<GeoLib::Point*> const& stations =
            *dynamic_cast<BaseItem*>(lists[i]->getItem())->getStations();

        for (std::size_t i = 0; i < stations.size(); i++)
        {
            auto* stratView = new StratWindow(
                static_cast<GeoLib::StationBorehole*>(stations[i]),
                &colorLookupTable);
            stratView->setAttribute(Qt::WA_DeleteOnClose);
            stratView->show();
            stratView->stationView->saveAsImage(QString::fromStdString(
                static_cast<GeoLib::StationBorehole*>(stations[i])->getName()) + ".jpg");
            stratView->close();
        }
    }
}
