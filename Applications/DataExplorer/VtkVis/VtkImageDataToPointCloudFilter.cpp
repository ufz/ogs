/**
 * \file
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkImageDataToPointCloudFilter.h"

#include <algorithm>
#include <vector>

#include <vtkIdList.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkLine.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

vtkStandardNewMacro(VtkImageDataToPointCloudFilter);

VtkImageDataToPointCloudFilter::VtkImageDataToPointCloudFilter()
    : Gamma(1.0),
      PointScaleFactor(1.0),
      MinHeight(0),
      MaxHeight(1000),
      MinNumberOfPointsPerCell(1),
      MaxNumberOfPointsPerCell(20),
      IsLinear(true)
{
}

void VtkImageDataToPointCloudFilter::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}

int VtkImageDataToPointCloudFilter::FillInputPortInformation(int, vtkInformation* info)
{
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    return 1;
}

int VtkImageDataToPointCloudFilter::RequestData(
    vtkInformation*,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector)
{
    vtkDebugMacro(<< "Executing VtkImageDataToPointCloudFilter");

    vtkInformation* input_info = inputVector[0]->GetInformationObject(0);
    vtkInformation* output_info = outputVector->GetInformationObject(0);
    vtkImageData* input = vtkImageData::SafeDownCast(
        input_info->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData* output = vtkPolyData::SafeDownCast(
        output_info->Get(vtkDataObject::DATA_OBJECT()));

    void* pixvals = input->GetScalarPointer();
    int const n_comp = input->GetNumberOfScalarComponents();
    if (n_comp < 1)
    {
        vtkDebugMacro("Error reading number of components.");
    }
    if (n_comp > 2)
    {
        vtkDebugMacro(
            "RGB colours detected. Using only first channel for computation.");
    }

    vtkIdType const n_points = input->GetNumberOfPoints();
    if (n_points == 0)
    {
        vtkDebugMacro("No data found!");
        return -1;
    }

    double spacing[3];
    input->GetSpacing(spacing);
    double range[2];
    vtkPointData* input_data = input->GetPointData();
    input_data->GetArray(0)->GetRange(range);

    std::vector<std::size_t> density;
    density.reserve(n_points);
    for (vtkIdType i = 0; i < n_points; ++i)
    {
        // Skip transparent pixels
        if (n_comp == 2 || n_comp == 4)
        {
            if (((float*)pixvals)[(i + 1) * n_comp - 1] < 0.00000001f)
            {
                density.push_back(0);
                continue;
            }
        }
        float const val(((float*)pixvals)[i * n_comp]);
        double const calc_gamma = (IsLinear) ? 1 : Gamma;
        double const pnts_per_cell = interpolate(range[0], range[1], val, calc_gamma);
        density.push_back(static_cast<std::size_t>(
            std::floor(pnts_per_cell * GetPointScaleFactor())));
    }

    // Allocate the space needed for output objects
    std::size_t const sum = std::accumulate(density.begin(), density.end(), 0);
    vtkSmartPointer<vtkPoints> new_points = vtkSmartPointer<vtkPoints>::New();
    new_points->SetNumberOfPoints(sum);
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    cells->Allocate(sum);

    double const half_cellsize(spacing[0] / 2.0);
    std::size_t pnt_idx(0);
    for (std::size_t i = 0; i < static_cast<std::size_t>(n_points); ++i)
    {
        if (density[i] == 0)
            continue;
        double p[3];
        input->GetPoint(i, p);

        MathLib::Point3d min_pnt{std::array<double, 3>{
            {p[0] - half_cellsize, p[1] - half_cellsize, MinHeight}}};
        MathLib::Point3d max_pnt{std::array<double, 3>{
            {p[0] + half_cellsize, p[1] + half_cellsize, MaxHeight}}};
        createPoints(new_points, cells, pnt_idx, density[i], min_pnt, max_pnt);
        pnt_idx += density[i];
    }

    output->SetPoints(new_points);
    output->SetVerts(cells);
    output->Squeeze();

    vtkDebugMacro(<< "Created " << new_points->GetNumberOfPoints() << " points.");
    return 1;
}

void VtkImageDataToPointCloudFilter::createPoints(
    vtkSmartPointer<vtkPoints>& points,
    vtkSmartPointer<vtkCellArray>& cells,
    std::size_t pnt_idx,
    std::size_t n_points,
    MathLib::Point3d const& min_pnt,
    MathLib::Point3d const& max_pnt) const
{
    for (std::size_t i = 0; i < n_points; ++i)
    {
        double p[3];
        p[0] = getRandomNumber(min_pnt[0], max_pnt[0]);
        p[1] = getRandomNumber(min_pnt[1], max_pnt[1]);
        p[2] = getRandomNumber(min_pnt[2], max_pnt[2]);
        points->SetPoint(pnt_idx + i, p);
        cells->InsertNextCell(1);
        cells->InsertCellPoint(pnt_idx + i);
    }
}

double VtkImageDataToPointCloudFilter::getRandomNumber(double const& min, double const& max) const
{
    return ((double)std::rand() / RAND_MAX) * (max - min) + min;
}

std::size_t VtkImageDataToPointCloudFilter::interpolate(double low,
                                                        double high,
                                                        double p,
                                                        double gamma) const
{
    assert(p >= low && p <= high);
    double const r = (p - low) / (high - low);
    return static_cast<std::size_t>(
        (MaxNumberOfPointsPerCell - MinNumberOfPointsPerCell) * std::pow(r, gamma) +
         MinNumberOfPointsPerCell);
}
