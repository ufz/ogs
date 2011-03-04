/**
 * \file VtkSelectionFilter.cpp
 * 2011/02/09 KR Initial implementation
 *
 */

// ** VTK INCLUDES **
#include "VtkSelectionFilter.h"

#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>

vtkStandardNewMacro(VtkSelectionFilter);
vtkCxxRevisionMacro(VtkSelectionFilter, "$Revision: 6995 $");


VtkSelectionFilter::VtkSelectionFilter()
: _threshold(0.0), _ifSmaller(0.0)
{
}

VtkSelectionFilter::~VtkSelectionFilter()
{
}

void VtkSelectionFilter::PrintSelf( ostream& os, vtkIndent indent )
{
	this->Superclass::PrintSelf(os,indent);
	os << indent << "== VtkSelectionFilter ==" << endl;
}

int VtkSelectionFilter::RequestData( vtkInformation*, 
							         vtkInformationVector** inputVector, 
								     vtkInformationVector* outputVector )
{
	if (this->_selection.empty()) 
	{
		std::cout << "VtkSelectionFilter - Error: Selection array is empty..." << std::endl;
		return 0;
	}
	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkUnstructuredGrid *input = vtkUnstructuredGrid::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	
	vtkSmartPointer<vtkDoubleArray> colors = vtkSmartPointer<vtkDoubleArray>::New();
		colors->SetNumberOfComponents(1);
		colors->SetName("Selection");

	size_t nCells = input->GetNumberOfCells();
	size_t arrayLength = this->_selection.size();
	if (nCells>arrayLength)
		std::cout << "VtkSelectionFilter - Warning: Number of cells exceeds selection array length. Surplus cells won't be examined." << std::endl;

	for (size_t i=0; i<arrayLength; i++)
	{
		colors->InsertNextValue(_selection[i]);
	}

	vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
		output->CopyStructure(input);
		output->GetPointData()->PassData(input->GetPointData());
		output->GetCellData()->PassData(input->GetCellData());
		output->GetCellData()->AddArray(colors);
		output->GetCellData()->SetActiveScalars("Selection");

	return 1;
}

void VtkSelectionFilter::SetSelectionArray(std::vector<double> selection, double threshold, bool ifSmaller)
{
	this->_selection = selection;
	this->_threshold = threshold;
	this->_ifSmaller = ifSmaller;
}