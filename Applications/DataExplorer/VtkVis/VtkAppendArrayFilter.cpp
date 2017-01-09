/**
 * \file
 * \author Karsten Rink
 * \date   2011-02-09
 * \brief  Implementation of the VtkAppendArrayFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** VTK INCLUDES **
#include "VtkAppendArrayFilter.h"

#include <logog/include/logog.hpp>

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
        ERR("VtkAppendArrayFilter::RequestData(): Selection array is empty.");
        return 0;
    }
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkUnstructuredGrid* input =
            vtkUnstructuredGrid::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkSmartPointer<vtkDoubleArray> colors = vtkSmartPointer<vtkDoubleArray>::New();
    colors->SetNumberOfComponents(1);
    std::size_t arrayLength = this->_array.size();
    colors->SetNumberOfValues(arrayLength);
    colors->SetName("Selection");

    std::size_t nCells = input->GetNumberOfCells();
    if (nCells > arrayLength)
        WARN("VtkAppendArrayFilter::RequestData(): Number of cells exceeds selection array length. Surplus cells won't be examined.");

    for (std::size_t i = 0; i < arrayLength; i++)
        colors->SetValue(i, _array[i]);

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
