/**
 * \file
 * \author Lars Bilke
 * \date   2010-10-06
 * \brief  Implementation of the VtkImageDataToPolyDataFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkImageDataToLinePolyDataFilter.h"

#include <vtkIdList.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkLine.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

vtkStandardNewMacro(VtkImageDataToLinePolyDataFilter);

VtkImageDataToLinePolyDataFilter::VtkImageDataToLinePolyDataFilter() :
        ImageSpacing(0.0)
{
    this->SetLengthScaleFactor(1.0);
}

VtkImageDataToLinePolyDataFilter::~VtkImageDataToLinePolyDataFilter() = default;

void VtkImageDataToLinePolyDataFilter::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
    os << indent << "LengthScaleFactor: " << this->LengthScaleFactor << "\n";
}

int VtkImageDataToLinePolyDataFilter::FillInputPortInformation(int, vtkInformation* info)
{
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    return 1;
}

int VtkImageDataToLinePolyDataFilter::RequestData(vtkInformation*,
                                                  vtkInformationVector** inputVector,
                                                  vtkInformationVector* outputVector)
{
    vtkDebugMacro(<< "Executing VtkImageDataToPolyDataFilter");

    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkImageData* input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData* output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    void* inScalarPtr = input->GetScalarPointer();
    int numScalarComponents = input->GetNumberOfScalarComponents();

    double spacing[3];
    input->GetSpacing(spacing);
    this->ImageSpacing = spacing[0];

    int dimensions[3];
    input->GetDimensions(dimensions);

    // Skip execution if there are no points
    vtkIdType numPts = input->GetNumberOfPoints();
    if (numPts < 1)
    {
        vtkDebugMacro("No data to extract lines!");
        return 1;
    }

    // Allocate the space needed for the output cells.
    output->Allocate(numPts);

    // Allocate space for a new set of points
    vtkSmartPointer<vtkPoints> newPts = vtkSmartPointer<vtkPoints>::New();
    newPts->SetNumberOfPoints(numPts * 2);

    // Allocate space for the data associated with the new set of points
    vtkPointData* inPD = input->GetPointData();
    vtkPointData* outPD = output->GetPointData();
    outPD->CopyAllocate(inPD, numPts * 16, numPts);

    // Compute scaling factor that max height is 0.1 * longest image dimension
    double range[2];
    inPD->GetArray(0)->GetRange(range);
    float scalingFactor = (std::max(dimensions[0], dimensions[1]) * spacing[0] * 0.1)
                          / std::max(range[0], range[1]);

    double dir[3] = {0, 0, 1};

    // Traverse all points creating another point with scalar distance in Z direction
    for (vtkIdType ptId = 0; ptId < numPts; ++ptId)
    {
        // Skip translucent pixels
        float opacity = ((float*)inScalarPtr)[ptId * numScalarComponents + 1];
        if (opacity < 0.00000001f)
            continue;

        // Compute length of the new line (scalar * LengthScaleFactor)
        float length = ((float*)inScalarPtr)[ptId * numScalarComponents]
                       * scalingFactor * this->LengthScaleFactor;
        //float length = (((unsigned char*)inScalarPtr)[ptId * numScalarComponents]* this->LengthScaleFactor > 50000) ? 50000 : ((unsigned char*)inScalarPtr)[ptId * numScalarComponents]* this->LengthScaleFactor;

        // Skip this line if length is zero
        if (length < 0.00000001f)
            continue;

        // Get the old point location
        double p[3];
        input->GetPoint(ptId, p);

        // Compute the new point location
        double newPt[3];
        for(std::size_t i = 0; i < 3; ++i)
            newPt[i] = p[i] + dir[i] * length;

        // Copy the old point
        newPts->SetPoint(ptId * 2, p);
        outPD->CopyData(inPD, ptId, ptId * 2);

        // Create the new point
        newPts->SetPoint(ptId * 2 + 1, newPt);
        outPD->CopyData(inPD, ptId, ptId * 2 + 1);

        // Create the line
        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, ptId * 2);
        line->GetPointIds()->SetId(1, ptId * 2 + 1);
        output->InsertNextCell(line->GetCellType(), line->GetPointIds());
    }

    // Store the new set of points in the output
    output->SetPoints(newPts);
    output->GetPointData()->GetArray(0)->SetName("Colors");

    // Avoid keeping extra memory around
    output->Squeeze();

    vtkDebugMacro(<< "Created: "
                  << newPts->GetNumberOfPoints() << " points, "
                  << output->GetNumberOfCells() << " lines");

    return 1;
}
