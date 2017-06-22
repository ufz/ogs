/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Definition of the StationTreeView class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <QContextMenuEvent>
#include <QTreeView>
#include "GeoLib/Station.h"

class vtkPolyDataAlgorithm;

/**
 * \brief A view for the StationTreeModel
 * \sa StationTreeModel, ModelTreeItem
 */
class StationTreeView : public QTreeView
{
    Q_OBJECT

public:
    /// Constructor
    StationTreeView(QWidget* parent = nullptr);

    /// Update the view to visualise changes made to the underlying data
    void updateView();

public slots:
    //void filterStations(std::vector<PropertyBounds> bounds);

protected slots:
    /// Instructions if the selection of items in the view has changed.
    void selectionChanged(const QItemSelection& selected,
                          const QItemSelection& deselected) override;

    /// Instructions if the selection of items in the view has changed by events outside the view (i.e. by actions made in the visualisation).
    void selectionChangedFromOutside(const QItemSelection &selected,
                                     const QItemSelection &deselected);

private:
    /// Actions to be taken after a right mouse click is performed in the station view.
    void contextMenuEvent(QContextMenuEvent* e) override;

    /// Create image files from all stratigraphies in a borehole vector
    void writeStratigraphiesAsImages(QString listName);

private slots:
    void addStationList();
    void displayStratigraphy();
    void exportList();
    void exportStation();
    void removeStationList();
    /// Calls a SetNameDialog.
    void setNameForElement();
    void writeToFile();
    void showDiagramPrefsDialog();

signals:
    void enableSaveButton(bool);
    void enableRemoveButton(bool);
    void geoItemSelected(const vtkPolyDataAlgorithm*, int);
    void removeGeoItemSelection();
    void itemSelectionChanged(const QItemSelection & selected,
                              const QItemSelection & deselected);
    void openStationListFile(int);
    void propertiesDialogRequested(std::string name);
    void requestNameChangeDialog(const std::string&, std::size_t);
    void stationListExportRequested(std::string listName, std::string fileName);
    void stationListRemoved(std::string name);
    void stationListSaved(QString listName, QString fileName);
    void diagramRequested(QModelIndex&);
};
