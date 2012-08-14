/**
 * \file StationTreeModel.h
 * KR Initial implementation
 */

#ifndef QSTATIONTREEMODEL_H
#define QSTATIONTREEMODEL_H

#include <vector>

#include "ModelTreeItem.h"
#include "Point.h"
#include "TreeModel.h"
#include <vtkPolyDataAlgorithm.h>

namespace GEOLIB
{
class Station;
class StationBorehole;
}

class QString;
class QModelIndex;
class PropertyBounds;

/**
 * \brief A model for the StationTreeView implementing a tree as a double-linked list.
 *
 * A model for the StationTreeView implementing a tree as a double-linked list.
 * In addition to a simple TreeModel each item also contains a 2D / 3D GraphicsItem for visualization.
 * \sa TreeModel, StationTreeView, TreeItem, ModelTreeItem
 */
class StationTreeModel : public TreeModel
{
	Q_OBJECT

public:
	StationTreeModel(QObject* parent = 0);
	~StationTreeModel();

	void addStationList(QString listName, const std::vector<GEOLIB::Point*>* stations);
	void filterStations(const std::string &name,
	                    const std::vector<GEOLIB::Point*>* stations,
	                    const std::vector<PropertyBounds> &bounds);
	const std::vector<ModelTreeItem*> &getLists() { return _lists; }
	QModelIndex index(int row, int column, const QModelIndex &parent = QModelIndex()) const;
	//BaseItem* itemFromIndex( const QModelIndex& index ) const;
	void removeStationList(QModelIndex index);
	void removeStationList(const std::string &name);
	GEOLIB::Station* stationFromIndex( const QModelIndex& index, QString &listName ) const;
	vtkPolyDataAlgorithm* vtkSource(const std::string &name) const;

private:
	std::vector<ModelTreeItem*> _lists;
};

#endif //QSTATIONTREEMODEL_H
