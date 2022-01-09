/**
 * \file
 * \author Karsten Rink
 * \date   2010-04-21
 * \brief  Implementation of the VtkColorByHeightFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** VTK INCLUDES **
#include "VtkColorByHeightFilter.h"

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkLookupTable.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>

#include "BaseLib/Logging.h"
#include "VtkColorLookupTable.h"

vtkStandardNewMacro(VtkColorByHeightFilter);

VtkColorByHeightFilter::VtkColorByHeightFilter()
    : ColorLookupTable(VtkColorLookupTable::New()), _tableRange{0.0, 0.0}
{
}

VtkColorByHeightFilter::~VtkColorByHeightFilter() = default;

void VtkColorByHeightFilter::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);

    double range[2];
    ColorLookupTable->GetTableRange(range);
    os << indent << "== VtkColorByHeightFilter ==" << endl;
    os << indent << "Range: " << range[0] << "-" << range[1] << endl;
    os << indent << "Interpolation Type:"
       << static_cast<int>(ColorLookupTable->getInterpolationType()) << endl;
}

vtkMTimeType VtkColorByHeightFilter::GetMTime()
{
    vtkMTimeType t1 = this->Superclass::GetMTime();
    if (this->ColorLookupTable)
    {
        vtkMTimeType const t2 = this->ColorLookupTable->GetMTime();
        return std::max(t1, t2);
    }
    return t1;
}
int VtkColorByHeightFilter::RequestData(vtkInformation* /*request*/,
                                        vtkInformationVector** inputVector,
                                        vtkInformationVector* outputVector)
{
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkPolyData* input =
        vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkSmartPointer<vtkDoubleArray> colors =
        vtkSmartPointer<vtkDoubleArray>::New();
    colors->SetNumberOfComponents(1);
    std::size_t nPoints = input->GetNumberOfPoints();
    colors->SetNumberOfValues(nPoints);
    colors->SetName("Colors");

    // Inserts height values as a new scalar array
    for (std::size_t i = 0; i < nPoints; i++)
    {
        double p[3];
        input->GetPoint(i, p);
        colors->SetValue(i, p[2]);
    }

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkPolyData* output =
        vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    output->CopyStructure(input);
    output->GetPointData()->PassData(input->GetPointData());
    output->GetCellData()->PassData(input->GetCellData());
    output->GetPointData()->AddArray(colors);
    output->GetPointData()->SetActiveScalars("Colors");

    return 1;
}

void VtkColorByHeightFilter::SetTableRange(double min, double max)
{
    if (min < max)
    {
        this->_tableRange[0] = min;
        this->_tableRange[1] = max;
        this->ColorLookupTable->SetTableRange(min, max);
    }
    else
    {
        ERR("VtkColorByHeightFilter::SetLimits(min, max) - Limits not changed "
            "because min value > max value.");
    }
}

void VtkColorByHeightFilter::SetTableRangeScaling(double scale)
{
    this->_tableRangeScaling = scale;
    this->ColorLookupTable->SetTableRange(this->_tableRange[0] * scale,
                                          this->_tableRange[1] * scale);
}
