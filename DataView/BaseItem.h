/**
 * \file BaseItem.h
 * 20/01/2010 KR Initial implementation
 *
 */
#ifndef BASEITEM_H
#define BASEITEM_H

#include "Point.h"

#include "VtkStationSource.h"
#include <QModelIndex>
#include <vtkPolyDataAlgorithm.h>

/**
 * \brief A BaseItem contains additional Information about a subtree in the StationTreeModel.
 *
 * It is used for list names in the StationTreeModel and it contains the
 * vtkObject for visualisation of the whole list in 3D.
 */
class BaseItem
{
public:
	BaseItem(const QString &listName, const std::vector<GEOLIB::Point*>* stations = NULL )
		: _stations(stations), _vtkSource(VtkStationSource::New())
	{
		// create the vtk-object for 3d-visualisation of this list
		static_cast<VtkStationSource*>(_vtkSource)->setStations(stations);
		static_cast<VtkStationSource*>(_vtkSource)->SetName(listName);
	}

	~BaseItem()
	{
		_vtkSource->Delete();
	}

	/// Returns the associated QModelIndex which belongs to a Qt model
	QModelIndex modelIndex() const { return _modelIndex; }

	/// Sets the model index
	void setModelIndex( QModelIndex index ) { _modelIndex = index; }

	const std::vector<GEOLIB::Point*>* getStations() { return _stations; }

	/// Returns the Vtk polydata source object
	vtkPolyDataAlgorithm* vtkSource() const { return _vtkSource; }

private:
	QModelIndex _modelIndex;
	const std::vector<GEOLIB::Point*>* _stations;

	/// The Vtk data source object. This is the starting point for a Vtk data
	/// visualization pipeline.
	vtkPolyDataAlgorithm* _vtkSource;
};

#endif //BASEITEM_H
