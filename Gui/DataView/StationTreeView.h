/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Definition of the StationTreeView class.
 *
 * \copyright
 * Copyright (c)  2013, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef QSTATIONTREEVIEW_H
#define QSTATIONTREEVIEW_H

#include "PropertyBounds.h"
#include <QContextMenuEvent>
#include <QTreeView>

#include "Station.h"

/**
 * \brief A view for the StationTreeModel with a number of properties adequate for this kind of data
 * \sa StationTreeModel, ModelTreeItem
 */
class StationTreeView : public QTreeView
{
	Q_OBJECT

public:
	/// Constructor
	StationTreeView(QWidget* parent = 0);

	/// Update the view to visualise changes made to the underlying data
	void updateView();

public slots:
	//void filterStations(std::vector<PropertyBounds> bounds);

protected slots:
	/// Instructions if the selection of items in the view has changed.
	void selectionChanged(const QItemSelection &selected, const QItemSelection &deselected);

	/// Instructions if the selection of items in the view has changed by events outside the view (i.e. by actions made in the visualisation).
	void selectionChangedFromOutside(const QItemSelection &selected,
	                                 const QItemSelection &deselected);

private:
	/// Actions to be taken after a right mouse click is performed in the station view.
	void contextMenuEvent( QContextMenuEvent* e );

	/// Create image files from all stratigraphies in a borehole vector
	void writeStratigraphiesAsImages(QString listName);

private slots:
	void addStationList();
	void displayStratigraphy();
	void exportList();
	void exportStation();
	void removeStationList();
	void writeToFile();
	void showPropertiesDialog();
	void showDiagramPrefsDialog();

signals:
	void enableSaveButton(bool);
	void enableRemoveButton(bool);
	void itemSelectionChanged(const QItemSelection & selected,
	                          const QItemSelection & deselected);
	void openStationListFile(int);
	void propertiesDialogRequested(std::string name);
	void stationListExportRequested(std::string listName, std::string fileName);
	void stationListRemoved(std::string name);
	void stationListSaved(QString listName, QString fileName);
	void diagramRequested(QModelIndex&);
};

#endif //QSTATIONTREEVIEW_H
