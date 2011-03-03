/**
 * \file GeoTreeView.h
 * 2011/02/07 KR Initial implementation
 */

#ifndef GEOTREEVIEW_H
#define GEOTREEVIEW_H

#include <QTreeView>
#include <QContextMenuEvent>
#include "GeoType.h"

/**
 * \brief A view for the GeoTreeModel 
 * \sa GeoTreeModel, GeoTreeItem
 */
class GeoTreeView : public QTreeView
{
	Q_OBJECT

public:
	/// Constructor
	GeoTreeView(QWidget* parent = 0);

	/// Update the view to visualise changes made to the underlying data
	void updateView();

protected slots:
	/// Instructions if the selection of items in the view has changed.
	void selectionChanged(const QItemSelection &selected, const QItemSelection &deselected);

	/// Instructions if the selection of items in the view has changed by events outside the view (i.e. by actions made in the visualisation).
	void selectionChangedFromOutside(const QItemSelection &selected, const QItemSelection &deselected);

private:
	/// Actions to be taken after a right mouse click is performed in the station view.
	void contextMenuEvent( QContextMenuEvent* e );

private slots:
	void on_Clicked(QModelIndex idx);
	/// Allows to add FEM Conditions to add to Geometry
	void addFEMConditions();
	/// Calls a LineEditDialog.
	void connectPolylines();
	/// Saves a geometry in a file.
	void writeToFile() const;
	/// Removes a whole geometry or parts of it.
	void removeList();
	
signals:
	void itemSelectionChanged(const QItemSelection & selected, const QItemSelection & deselected);
	void listRemoved(std::string name, GEOLIB::GEOTYPE);
	void loadFEMCondFileRequested(std::string);
	void saveToFileRequested(QString, QString) const;
	void requestLineEditDialog(const std::string&);
};

#endif //GEOTREEVIEW_H
