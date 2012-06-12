/**
 * \file VtkCompositeThresholdFilter.cpp
 * 25/10/2010 LB Initial implementation
 *
 * Implementation of VtkCompositeThresholdFilter class
 */

// ** INCLUDES **
#include "VtkCompositeThresholdFilter.h"

#include <vtkCellData.h>
#include <vtkThreshold.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIntArray.h>
#include <vtkSmartPointer.h>

#include <limits>

VtkCompositeThresholdFilter::VtkCompositeThresholdFilter( vtkAlgorithm* inputAlgorithm )
	: VtkCompositeFilter(inputAlgorithm)
{
	this->init();
}

VtkCompositeThresholdFilter::~VtkCompositeThresholdFilter()
{
}

void VtkCompositeThresholdFilter::init()
{
	// Set meta information about input and output data types
	this->_inputDataObjectType = VTK_DATA_SET;
	this->_outputDataObjectType = VTK_UNSTRUCTURED_GRID;

	// Because this is the only filter here we cannot use vtkSmartPointer
	vtkThreshold* threshold = vtkThreshold::New();
	threshold->SetInputConnection(_inputAlgorithm->GetOutputPort());

	// Note: There is no need to select the input array because it is
	//       automatically selected.

	// Sets a filter property which will be user editable
	threshold->SetSelectedComponent(0);

	// Setting the threshold to min / max values to ensure that the whole data
	// is first processed. This is needed for correct lookup table generation.
	const double dMin = std::numeric_limits<double>::min();
	const double dMax = std::numeric_limits<double>::max();
	threshold->ThresholdBetween(dMin, dMax);

	// Create a list for the ThresholdBetween (vector) property.
	QList<QVariant> thresholdRangeList;
	// Insert the values (same values as above)
	thresholdRangeList.push_back(dMin);
	thresholdRangeList.push_back(dMax);
	// Put that list in the property map
	(*_algorithmUserVectorProperties)["Threshold Between"] = thresholdRangeList;

	// Make a new entry in the property map for the SelectedComponent property
	(*_algorithmUserProperties)["Selected Component"] = 0;

	// Must all scalars match the criterium
	threshold->SetAllScalars(1);
	(*_algorithmUserProperties)["Evaluate all points"] = true;

	// The threshold filter is last one and so it is also the _outputAlgorithm
	_outputAlgorithm = threshold;
}

void VtkCompositeThresholdFilter::SetUserProperty( QString name, QVariant value )
{
	VtkAlgorithmProperties::SetUserProperty(name, value);

	// Use the same name as in init()
	if (name.compare("Selected Component") == 0)
		// Set the property on the algorithm
		static_cast<vtkThreshold*>(_outputAlgorithm)->SetSelectedComponent(value.toInt());
	else if (name.compare("Evaluate all points") == 0)
		static_cast<vtkThreshold*>(_outputAlgorithm)->SetAllScalars(value.toBool());
}

void VtkCompositeThresholdFilter::SetUserVectorProperty( QString name, QList<QVariant> values )
{
	VtkAlgorithmProperties::SetUserVectorProperty(name, values);

	// Use the same name as in init()
	if (name.compare("Threshold Between") == 0)
		// Set the vector property on the algorithm
		static_cast<vtkThreshold*>(_outputAlgorithm)->ThresholdBetween(
		        values[0].toDouble(), values[1].toDouble());
}
