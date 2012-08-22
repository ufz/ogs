/**
 * \file GeoTreeView.h
 * 2011/02/07 KR Initial implementation
 */

#ifndef GEOTREEVIEW_H
#define GEOTREEVIEW_H

#include "GeoType.h"
#include <QContextMenuEvent>
#include <QTreeView>

class vtkPolyDataAlgorithm;

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
	void selectionChangedFromOutside(const QItemSelection &selected,
	                                 const QItemSelection &deselected);

private:
	/// Actions to be taken after a right mouse click is performed in the station view.
	void contextMenuEvent( QContextMenuEvent* e );
	/// Calls a FEMConditionSetupDialog.
	void setElementAsCondition(bool set_on_points = false);

private slots:
	/// Allows to add FEM Conditions to a process
	void loadFEMConditions();
	void on_Clicked(QModelIndex idx);
	/// Calls a LineEditDialog.
	void connectPolylines();
	/// Calls a SetNameDialog.
	void setNameForElement();
	void setObjectAsCondition() { setElementAsCondition(false); };
	void setObjectPointsAsCondition() { setElementAsCondition(true); };
	/// Saves a geometry in a file.
	void writeToFile() const;
	/// Removes a whole geometry or parts of it.
	void removeList();
	/// Saves FEM Conditions associated with the given geometry
	//void saveFEMConditions();

signals:
	void geoItemSelected(const vtkPolyDataAlgorithm*, int);
	void removeGeoItemSelection();
	//void itemSelectionChanged(const QItemSelection & selected, const QItemSelection & deselected);
	void listRemoved(std::string name, GeoLib::GEOTYPE);
	void loadFEMCondFileRequested(std::string);
	void saveToFileRequested(QString, QString) const;
	void requestCondSetupDialog(const std::string&, const GeoLib::GEOTYPE, const size_t, bool on_points);
	void requestLineEditDialog(const std::string&);
	void requestNameChangeDialog(const std::string&, const GeoLib::GEOTYPE, const size_t);
	//void saveFEMConditionsRequested(QString, QString);
};

#endif //GEOTREEVIEW_H
