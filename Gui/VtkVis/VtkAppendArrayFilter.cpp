/**
 * \file
 * \author Karsten Rink
 * \date   2011-02-09
 * \brief  Implementation of the VtkAppendArrayFilter class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** VTK INCLUDES **
#include "VtkAppendArrayFilter.h"

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(VtkAppendArrayFilter);
vtkCxxRevisionMacro(VtkAppendArrayFilter, "$Revision$");

VtkAppendArrayFilter::VtkAppendArrayFilter()
{
}

VtkAppendArrayFilter::~VtkAppendArrayFilter()
{
}

void VtkAppendArrayFilter::PrintSelf( ostream& os, vtkIndent indent )
{
	this->Superclass::PrintSelf(os,indent);
	os << indent << "== VtkAppendArrayFilter ==" << endl;
}

int VtkAppendArrayFilter::RequestData( vtkInformation*,
                                     vtkInformationVector** inputVector,
                                     vtkInformationVector* outputVector )
{
	if (this->_array.empty())
	{
		std::cout << "VtkAppendArrayFilter - Error: Selection array is empty..." << std::endl;
		return 0;
	}
	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
	vtkUnstructuredGrid* input =
	        vtkUnstructuredGrid::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkSmartPointer<vtkDoubleArray> colors = vtkSmartPointer<vtkDoubleArray>::New();
	colors->SetNumberOfComponents(1);
	colors->SetName("Selection");

	size_t nCells = input->GetNumberOfCells();
	size_t arrayLength = this->_array.size();
	if (nCells > arrayLength)
		std::cout <<
		"VtkAppendArrayFilter - Warning: Number of cells exceeds selection array length. Surplus cells won't be examined."
		          << std::endl;

	for (size_t i = 0; i < arrayLength; i++)
		colors->InsertNextValue(_array[i]);

	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	vtkUnstructuredGrid* output =
	        vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
	output->CopyStructure(input);
	output->GetPointData()->PassData(input->GetPointData());
	output->GetCellData()->PassData(input->GetCellData());
	output->GetCellData()->AddArray(colors);
	output->GetCellData()->SetActiveScalars(_array_name.c_str());

	return 1;
}

void VtkAppendArrayFilter::SetArray(const std::string &array_name,
									const std::vector<double> &new_array)
{
	this->_array_name = array_name;
	this->_array = new_array;
}
