/**
 * \file VtkCompositeContourFilter.cpp
 * 2011/08/05 KR Initial implementation
 */

// ** INCLUDES **
#include "VtkCompositeContourFilter.h"

#include <vtkCellData.h>
#include <vtkContourFilter.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

VtkCompositeContourFilter::VtkCompositeContourFilter( vtkAlgorithm* inputAlgorithm )
	: VtkCompositeFilter(inputAlgorithm)
{
	this->init();
}

VtkCompositeContourFilter::~VtkCompositeContourFilter()
{
}

void VtkCompositeContourFilter::init()
{
	// Set meta information about input and output data types
	this->_inputDataObjectType = VTK_UNSTRUCTURED_GRID; //VTK_DATA_SET;
	this->_outputDataObjectType = VTK_UNSTRUCTURED_GRID;

	// Because this is the only filter here we cannot use vtkSmartPointer
	vtkContourFilter* contour = vtkContourFilter::New();
	contour->SetInputConnection(_inputAlgorithm->GetOutputPort());

	// Sets a filter vector property which will be user editable
	contour->GenerateValues(10, 0, 100);

	// Create a list for the ThresholdBetween (vector) property.
	QList<QVariant> contourRangeList;
	// Insert the values (same values as above)
	contourRangeList.push_back(0);
	contourRangeList.push_back(100);
	// Put that list in the property map
	(*_algorithmUserVectorProperties)["Range"] = contourRangeList;

	// Make a new entry in the property map for the "Number of Values" property
	(*_algorithmUserProperties)["Number of Contours"] = 10;

	// The threshold filter is last one and so it is also the _outputAlgorithm
	_outputAlgorithm = contour;
}

void VtkCompositeContourFilter::SetUserProperty( QString name, QVariant value )
{
	VtkAlgorithmProperties::SetUserProperty(name, value);

	// Use the same name as in init()
	if (name.compare("Number of Contours") == 0)
		static_cast<vtkContourFilter*>(_outputAlgorithm)->SetNumberOfContours(value.toInt());
}

void VtkCompositeContourFilter::SetUserVectorProperty( QString name, QList<QVariant> values )
{
	VtkAlgorithmProperties::SetUserVectorProperty(name, values);

	// Use the same name as in init()
	if (name.compare("Range") == 0)
		static_cast<vtkContourFilter*>(_outputAlgorithm)->GenerateValues(
		        VtkAlgorithmProperties::GetUserProperty("Number of Contours").toInt(),
		        values[0].toDouble(),
		        values[1].toDouble());
}
