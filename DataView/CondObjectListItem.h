/**
 * \file GeoObjectListItem.h
 * 2011/04/05 KR Initial implementation
 *
 */
#ifndef CONDOBJECTLISTITEM_H
#define CONDOBJECTLISTITEM_H


#include "FEMCondition.h"
#include "TreeItem.h"
#include <vtkPolyDataAlgorithm.h>
#include "VtkConditionSource.h"
#include <QModelIndex>

/**
 * \brief The CondObjectListItem is the TreeItem that contains the subtree for either initial conditions, 
 * boundary conditions source terms. This item also contains the vtk source-item for visualisation of this
 * information and the indices of the associated geometry-objects.
 * \sa TreeItem
 */
class CondObjectListItem : public TreeItem
{

public:
	/// Constructor for the TreeItem specifying FEM Conditions.
	CondObjectListItem(const QList<QVariant> &data, TreeItem *parent, FEMCondition::CondType type, const std::vector<GEOLIB::Point*> *points, const std::vector<GEOLIB::Polyline*> *lines, const std::vector<GEOLIB::Surface*> *surfaces)
		: TreeItem(data, parent), _vtkSource(VtkConditionSource::New()),  _type(type), _cond_vec(new std::vector<FEMCondition*>), _use_domain(new bool(false))
	{
		QString display_name = parent->data(0).toString().append(" - ").append(QString::fromStdString(FEMCondition::condTypeToString(type)));
		static_cast<VtkConditionSource*>(_vtkSource)->setData( points, _cond_vec, _use_domain);
		static_cast<VtkConditionSource*>(_vtkSource)->SetName( display_name );
	}

	~CondObjectListItem()
	{
		_vtkSource->Delete();
		delete _cond_vec;
		delete _use_domain;
	}

	/// Adds FEMCondtion for construction of VTK object.
	void addCondition(FEMCondition* cond) { 
		_cond_vec->push_back(cond); 
		_vtkSource->Modified();
	};


	/// Returns the type of geo-objects contained in the subtree of this item.
	FEMCondition::CondType getType() { return _type; };

	/// Returns the Vtk polydata source object
	vtkPolyDataAlgorithm* vtkSource() const { return _vtkSource; };

private:
	/// The Vtk data source object. This is the starting point for a Vtk data visualization pipeline.
	vtkPolyDataAlgorithm* _vtkSource;
	FEMCondition::CondType _type;
	std::vector<FEMCondition*> *_cond_vec;
	bool* _use_domain;
};

#endif //CONDOBJECTLISTITEM_H
