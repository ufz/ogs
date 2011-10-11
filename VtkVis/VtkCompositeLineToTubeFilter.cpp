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

	vtkTubeFilter* tubes = vtkTubeFilter::New();
	tubes->SetInputConnection(0, mergePoints->GetOutputPort(0));
	tubes->SetInputArrayToProcess(1,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"Stratigraphies");
	tubes->SetRadius(150);
	tubes->SetNumberOfSides(10);
	tubes->SetCapping(1);
	//tubes->SetVaryRadiusToVaryRadiusByScalar(); // KR radius changes with scalar

	(*_algorithmUserProperties)["Radius"] = 150.0;
	(*_algorithmUserProperties)["NumberOfSides"] = 6;
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
