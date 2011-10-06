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

	// TODO: Initalisation
	//vtkDataArray* activeScalar = vtkDataSet::SafeDownCast(threshold->GetInput())->GetCellData()->GetArray(0);
	//threshold->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS, activeScalar->GetName());
	//threshold->SetComponentModeToUseSelected();

//double* range = threshold->GetOutput()->GetScalarRange();
//std::cout << range[0] << ", " << range[1] << std::endl;

	// Sets a filter property which will be user editable
	threshold->SetSelectedComponent(0);

	// Sets a filter vector property which will be user editable
	threshold->ThresholdBetween(0, 100);

	// Create a list for the ThresholdBetween (vector) property.
	QList<QVariant> thresholdRangeList;
	// Insert the values (same values as above)
	thresholdRangeList.push_back(0);
	thresholdRangeList.push_back(100);
	// Put that list in the property map
	(*_algorithmUserVectorProperties)["Threshold Between"] = thresholdRangeList;

	// Make a new entry in the property map for the SelectedComponent property
	(*_algorithmUserProperties)["Selected Component"] = 0;

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
}

void VtkCompositeThresholdFilter::SetUserVectorProperty( QString name, QList<QVariant> values )
{
	VtkAlgorithmProperties::SetUserVectorProperty(name, values);

	// Use the same name as in init()
	if (name.compare("Threshold Between") == 0)
		//double* range = dynamic_cast<vtkUnstructuredGridAlgorithm*>(_outputAlgorithm)->GetOutput()->GetScalarRange();
		//std::cout << range[0] << ", " << range[1] << std::endl;
		// Set the vector property on the algorithm
		static_cast<vtkThreshold*>(_outputAlgorithm)->ThresholdBetween(
		        values[0].toInt(), values[1].toInt());
}
/*
   void VtkCompositeThresholdFilter::SetScalarRangeOnItem( double min, double max )
   {
    _item->SetScalarRange(min, max);
    emit requestViewUpdate();
   }

 */
