/**
 * \file VtkCompositeSelectionFilter.cpp
 * 2011/02/10 KR Initial implementation
 * 
 * Implementation of VtkCompositeSelectionFilter class
 */

// ** INCLUDES **
#include "VtkCompositeSelectionFilter.h"
#include "VtkSelectionFilter.h"

#include <vtkSmartPointer.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkThreshold.h>

VtkCompositeSelectionFilter::VtkCompositeSelectionFilter( vtkAlgorithm* inputAlgorithm )
: VtkCompositeFilter(inputAlgorithm)
{
	//this->init();
}

void VtkCompositeSelectionFilter::init()
{
	double thresholdValue(1.0);
	this->_inputDataObjectType = VTK_UNSTRUCTURED_GRID;
	this->_outputDataObjectType = VTK_UNSTRUCTURED_GRID;

	VtkSelectionFilter* selFilter = VtkSelectionFilter::New();
		selFilter->SetInputConnection(_inputAlgorithm->GetOutputPort());
		selFilter->SetSelectionArray(_selection, thresholdValue, true);
		selFilter->Update();

	vtkThreshold* threshold = vtkThreshold::New();
		threshold->SetInputConnection(selFilter->GetOutputPort());
		threshold->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS, "Selection");
		threshold->SetSelectedComponent(0);
		threshold->ThresholdByLower(thresholdValue);
		threshold->Update();
	(*_algorithmUserProperties)["Threshold"] = thresholdValue;

	_outputAlgorithm = threshold;
}

void VtkCompositeSelectionFilter::SetUserProperty( QString name, QVariant value )
{
	VtkAlgorithmProperties::SetUserProperty(name, value);

	if (name.compare("Threshold") == 0)
		// Set the vector property on the algorithm
		static_cast<vtkThreshold*>(_outputAlgorithm)->ThresholdByLower(value.toDouble());
}
