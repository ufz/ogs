/**
 * \file GeoObjectListItem.h
 * 2011/04/05 KR Initial implementation
 *
 */
#ifndef CONDOBJECTLISTITEM_H
#define CONDOBJECTLISTITEM_H

#include "FEMCondition.h"
#include "TreeItem.h"
#include "VtkConditionSource.h"
#include <QModelIndex>
#include <vtkPolyDataAlgorithm.h>
#include <vtkThresholdPoints.h>

/**
 * \brief The CondObjectListItem is the TreeItem that contains the subtree for either initial conditions,
 * boundary conditions source terms. This item also contains the vtk source-item for visualisation of this
 * information and the indices of the associated geometry-objects.
 * Upon creation the type of condition needs to be defined and the vector of points of the associated geometry
 * is needed for created of the vtk-object.
 * The actual FEM Condtions are added using the addCondition()-method.
 * \sa TreeItem
 */
class CondObjectListItem : public TreeItem
{
public:
	/// Constructor for the TreeItem specifying FEM Conditions.
	CondObjectListItem(const QList<QVariant> &data,
	                   TreeItem* parent,
	                   const FEMCondition::CondType type,
	                   const std::vector<GEOLIB::Point*>* points)
		: TreeItem(data, parent), _vtkSource(VtkConditionSource::New()),  _type(type),
		  _cond_vec(new std::vector<FEMCondition*>)
	{
		QString display_name = parent->data(0).toString().append(" - ").append(
		        QString::fromStdString(FEMCondition::condTypeToString(type)));
		static_cast<VtkConditionSource*>(_vtkSource)->setData( points, _cond_vec);
		static_cast<VtkConditionSource*>(_vtkSource)->SetName( display_name );
	}

	~CondObjectListItem()
	{
		_vtkSource->Delete();
		delete _cond_vec;
	}

	/// Adds FEMCondtion for construction of VTK object.
	void addCondition(FEMCondition* cond)
	{
		_cond_vec->push_back(cond);
		_vtkSource->Modified();
	}

	/// Returns the type of geo-objects contained in the subtree of this item.
	const FEMCondition::CondType getType() const { return _type; }

	/// Returns the Vtk polydata source object
	vtkPolyDataAlgorithm* vtkSource() const
	{
		return _vtkSource;
		/*
		   vtkThresholdPoints* threshold = vtkThresholdPoints::New();
		   threshold->SetInputConnection(_vtkSource->GetOutputPort());
		   threshold->ThresholdByUpper(-9998);
		   threshold->Update();
		   return threshold;
		 */
	}

private:
	/// The Vtk data source object. This is the starting point for a Vtk data visualization pipeline.
	vtkPolyDataAlgorithm* _vtkSource;
	FEMCondition::CondType _type;
	std::vector<FEMCondition*>* _cond_vec;
};

#endif //CONDOBJECTLISTITEM_H
