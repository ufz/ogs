/**
 * \file
 * \author Karsten Rink
 * \date   2011-02-07
 * \brief  Implementation of the GeoTreeView class.
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
#include <QSettings>

#include "GeoObjectListItem.h"
#include "GeoTreeItem.h"
#include "GeoTreeModel.h"
#include "GeoTreeView.h"
#include "OGSError.h"
#include "ImportFileTypes.h"
#include "LastSavedFileDirectory.h"


GeoTreeView::GeoTreeView(QWidget* parent) : QTreeView(parent)
{
}

void GeoTreeView::updateView()
{
    setAlternatingRowColors(true);
    //resizeColumnToContents(0);
    setColumnWidth(0,150);
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
        emit removeGeoItemSelection();

        const GeoObjectListItem* geo_object = dynamic_cast<GeoObjectListItem*>(tree_item->parentItem());
        if (geo_object) // geometry object
        {
            emit enableSaveButton(false);
            emit enableRemoveButton(false);
            emit geoItemSelected(geo_object->vtkSource(), tree_item->row());
        }
        else
        {
            if (!idx.parent().isValid()) // geometry item
            {
                emit enableSaveButton(true);
                emit enableRemoveButton(true);
            }
            else // line points or surface triangles
            {
                emit enableSaveButton(false);
                const auto* geo_type =
                    dynamic_cast<const GeoObjectListItem*>(tree_item);
                if (geo_type) // geometry list item
                    emit enableRemoveButton(true);
                else
                {
                    // highlight a point for an expanded polyline
                    auto* list_item = dynamic_cast<GeoObjectListItem*>(
                        tree_item->parentItem()->parentItem());
                    if (list_item && list_item->getType() == GeoLib::GEOTYPE::POLYLINE)
                        geoItemSelected(
                            dynamic_cast<GeoObjectListItem*>(tree_item->parentItem()->parentItem()->parentItem()->child(0))->vtkSource(),
                            tree_item->data(0).toInt());

                    // highlight a point for an expanded surface
                    list_item = dynamic_cast<GeoObjectListItem*>(tree_item->parentItem()->parentItem()->parentItem());
                    if (list_item && list_item->getType() == GeoLib::GEOTYPE::SURFACE)
                        geoItemSelected(
                            dynamic_cast<GeoObjectListItem*>(tree_item->parentItem()->parentItem()->parentItem()->parentItem()->child(0))->vtkSource(),
                            tree_item->data(0).toInt());
                    emit enableRemoveButton(false);
                }
            }

        }

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

    QTreeView::selectionChanged(selected, deselected);
}

void GeoTreeView::contextMenuEvent( QContextMenuEvent* event )
{
    QModelIndex index = this->selectionModel()->currentIndex();
    auto* item = static_cast<TreeItem*>(index.internalPointer());

    auto* list = dynamic_cast<GeoObjectListItem*>(item);
    QMenu menu;

    // The current index is a list of points/polylines/surfaces
    if (list != nullptr)
    {
        QAction* connectPlyAction(nullptr);
        if (list->getType() == GeoLib::GEOTYPE::POLYLINE)
        {
            connectPlyAction = menu.addAction("Connect Polylines...");
            connect(connectPlyAction, SIGNAL(triggered()), this,
                    SLOT(connectPolylines()));
        }
        menu.addSeparator();
        //QAction* removeAction = menu.addAction("Remove " + item->data(0).toString());
        //connect(removeAction, SIGNAL(triggered()), this, SLOT(removeList()));
    }
    else
    {
        if (!item)  // Otherwise sometimes it crashes when (unmotivated ;-) ) clicking in a treeview
            return;

        auto* parent = dynamic_cast<GeoObjectListItem*>(item->parentItem());

        // The current index refers to a geo-object
        if (parent != nullptr)
        {
            QMenu* cond_menu = new QMenu("Set FEM Condition");
            //menu.addMenu(cond_menu);
            QAction* addCondAction = cond_menu->addAction("On object...");
            QAction* addCondPointAction = cond_menu->addAction("On all points...");
            QAction* addNameAction = menu.addAction("Set name...");
            connect(addCondAction, SIGNAL(triggered()), this, SLOT(setObjectAsCondition()));
            connect(addNameAction, SIGNAL(triggered()), this, SLOT(setNameForElement()));

            if (parent->getType() == GeoLib::GEOTYPE::POINT)
                addCondPointAction->setEnabled(false);
            else
                connect(addCondPointAction, SIGNAL(triggered()), this, SLOT(setObjectPointsAsCondition()));

        }
        // The current index refers to the name of a geometry-object
        else if (item->childCount() > 0)
        {
            if (item->child(0)->data(0).toString().compare("Points") == 0) // clumsy way to find out
            {
                //QAction* saveAction = menu.addAction("Save geometry...");
                QAction* mapAction = menu.addAction("Map geometry...");
                //QAction* addCNDAction = menu.addAction("Load FEM Conditions...");
                //QAction* saveCondAction    = menu.addAction("Save FEM conditions...");
                menu.addSeparator();
                //QAction* removeAction = menu.addAction("Remove geometry");
                //connect(saveAction, SIGNAL(triggered()), this, SLOT(writeToFile()));
                connect(mapAction, SIGNAL(triggered()), this, SLOT(mapGeometry()));
                //connect(addCNDAction, SIGNAL(triggered()), this, SLOT(loadFEMConditions()));
                //connect(saveCondAction, SIGNAL(triggered()), this, SLOT(saveFEMConditions()));
                //connect(removeAction, SIGNAL(triggered()), this, SLOT(removeList()));
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

void GeoTreeView::addGeometry()
{
    emit openGeometryFile(ImportFileType::OGS_GEO);
}

void GeoTreeView::removeGeometry()
{
    QModelIndex index (this->selectionModel()->currentIndex());
    if (!index.isValid())
        OGSError::box("No geometry selected.");
    else
    {
        TreeItem* item = static_cast<GeoTreeModel*>(model())->getItem(index);
        auto* list = dynamic_cast<GeoObjectListItem*>(item);
        if (list) {
            emit listRemoved((item->parentItem()->data(
                                      0).toString()).toStdString(), list->getType());
        } else {
            emit listRemoved((item->data(0).toString()).toStdString(), GeoLib::GEOTYPE::SURFACE);
            emit listRemoved((item->data(0).toString()).toStdString(), GeoLib::GEOTYPE::POLYLINE);
            emit listRemoved((item->data(0).toString()).toStdString(), GeoLib::GEOTYPE::POINT);
        }

        if(this->selectionModel()->selectedIndexes().count() == 0)
        {
            emit enableSaveButton(false);
            emit enableRemoveButton(false);
        }
    }
}

void GeoTreeView::setElementAsCondition(bool set_on_points)
{
    const TreeItem* item = static_cast<GeoTreeModel*>(model())->getItem(
            this->selectionModel()->currentIndex());
    const std::size_t id = item->row();
    const GeoLib::GEOTYPE type = static_cast<GeoObjectListItem*>(item->parentItem())->getType();
    const std::string geometry_name = item->parentItem()->parentItem()->data(0).toString().toStdString();
    emit requestCondSetupDialog(geometry_name, type, id, set_on_points);
}

void GeoTreeView::setNameForElement()
{
    const TreeItem* item = static_cast<GeoTreeModel*>(model())->getItem(
            this->selectionModel()->currentIndex());
    const std::size_t id = item->row();
    const GeoLib::GEOTYPE type = static_cast<GeoObjectListItem*>(item->parentItem())->getType();
    const std::string geometry_name = item->parentItem()->parentItem()->data(0).toString().toStdString();
    emit requestNameChangeDialog(geometry_name, type, id);
}

void GeoTreeView::mapGeometry()
{
    TreeItem* item = static_cast<GeoTreeModel*>(model())->getItem(
            this->selectionModel()->currentIndex());
    std::string geo_name (item->data(0).toString().toStdString());
    emit geometryMappingRequested(geo_name);
}

void GeoTreeView::writeToFile() const
{
    QModelIndex index (this->selectionModel()->currentIndex());
    if (!index.isValid())
        OGSError::box("No geometry selected.");
    else
    {
        TreeItem* item = static_cast<GeoTreeModel*>(model())->getItem(index);
        QString file_type ("GeoSys geometry file (*.gml)");
#ifndef NDEBUG
         file_type.append(";;Legacy geometry file (*.gli)");
#endif // DEBUG
        QString geoName = item->data(0).toString();
        QString fileName = QFileDialog::getSaveFileName(
            nullptr, "Save geometry as",
            LastSavedFileDirectory::getDir() + geoName, file_type);
        if (!fileName.isEmpty())
        {
            LastSavedFileDirectory::setDir(fileName);
            emit saveToFileRequested(geoName, fileName);
        }
    }
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
    QString fileName =
        QFileDialog::getSaveFileName(nullptr, "Save FEM Conditions as", "",
                                     "OpenGeoSys FEM Condition file (*.cnd);; "
                                     "GeoSys Boundary Condition (*.bc);; "
                                     "GeoSys Initial Condition (*.ic);; "
                                     "GeoSys Source Condition (*.st)");
    emit saveFEMConditionsRequested(item->data(0).toString(), fileName);
}
*/
