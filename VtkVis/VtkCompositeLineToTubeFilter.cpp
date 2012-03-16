/**
 * \file VtkCompositeLineToTubeFilter.cpp
 * 18/11/2010 KR Initial implementation
 *
 * Implementation of VtkCompositeLineToTubeFilter class
 */

// ** INCLUDES **
#include "VtkCompositeLineToTubeFilter.h"

#include <vtkCleanPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTubeFilter.h>

VtkCompositeLineToTubeFilter::VtkCompositeLineToTubeFilter( vtkAlgorithm* inputAlgorithm )
	: VtkCompositeFilter(inputAlgorithm)
{
	this->init();
}

VtkCompositeLineToTubeFilter::~VtkCompositeLineToTubeFilter()
{
}

void VtkCompositeLineToTubeFilter::init()
{
	this->_inputDataObjectType = VTK_DATA_SET;
	this->_outputDataObjectType = VTK_POLY_DATA;

	// collapse coincident points
	vtkSmartPointer<vtkCleanPolyData> mergePoints = vtkSmartPointer<vtkCleanPolyData>::New();
	mergePoints->SetInputConnection(0, _inputAlgorithm->GetOutputPort(0));
	mergePoints->SetTolerance(0.0);
	mergePoints->ConvertLinesToPointsOn();

	double default_radius(GetInitialRadius());
	int default_number_of_sides(8);
	vtkTubeFilter* tubes = vtkTubeFilter::New();
	tubes->SetInputConnection(0, mergePoints->GetOutputPort(0));

	//tubes->SetInputArrayToProcess(1,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"StationValue");
	//tubes->SetVaryRadiusToVaryRadiusByScalar(); // KR radius changes with scalar

	tubes->SetInputArrayToProcess(1,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"Stratigraphies");
	tubes->SetRadius(default_radius);
	tubes->SetNumberOfSides(default_number_of_sides);
	tubes->SetCapping(1);

	(*_algorithmUserProperties)["Radius"] = default_radius;
	(*_algorithmUserProperties)["NumberOfSides"] = default_number_of_sides;
	(*_algorithmUserProperties)["Capping"] = true;

	_outputAlgorithm = tubes;
}

void VtkCompositeLineToTubeFilter::SetUserProperty( QString name, QVariant value )
{
	VtkAlgorithmProperties::SetUserProperty(name, value);

	if (name.compare("Radius") == 0)
		static_cast<vtkTubeFilter*>(_outputAlgorithm)->SetRadius(value.toDouble());
	else if (name.compare("NumberOfSides") == 0)
		static_cast<vtkTubeFilter*>(_outputAlgorithm)->SetNumberOfSides(value.toInt());
	else if (name.compare("Capping") == 0)
		static_cast<vtkTubeFilter*>(_outputAlgorithm)->SetCapping(value.toBool());
}

float VtkCompositeLineToTubeFilter::GetInitialRadius() const
{
	double bounding_box[6];
	static_cast<vtkPolyData*>(this->_inputAlgorithm->GetOutputDataObject(0))->GetBounds(bounding_box);
	double x_diff = fabs(bounding_box[0]-bounding_box[1]);
	double y_diff = fabs(bounding_box[2]-bounding_box[3]);
	double z_diff = fabs(bounding_box[4]-bounding_box[5]);

	double max = (x_diff > y_diff) ? x_diff : y_diff;
	max = (max > z_diff) ? max : z_diff;

	return max/200.0;
}
