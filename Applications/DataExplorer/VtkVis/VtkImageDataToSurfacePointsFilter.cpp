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
#include "VtkImageDataToSurfacePointsFilter.h"

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

vtkStandardNewMacro(VtkImageDataToSurfacePointsFilter);

VtkImageDataToSurfacePointsFilter::VtkImageDataToSurfacePointsFilter()
    : PointsPerPixel(20)
{
}

void VtkImageDataToSurfacePointsFilter::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}

int VtkImageDataToSurfacePointsFilter::FillInputPortInformation(int, vtkInformation* info)
{
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    return 1;
}

int VtkImageDataToSurfacePointsFilter::RequestData(
    vtkInformation*,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector)
{
    vtkDebugMacro(<< "Executing VtkImageDataToSurfacePointsFilter");

    vtkInformation* input_info = inputVector[0]->GetInformationObject(0);
    vtkInformation* output_info = outputVector->GetInformationObject(0);
    vtkImageData* input = vtkImageData::SafeDownCast(
        input_info->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData* output = vtkPolyData::SafeDownCast(
        output_info->Get(vtkDataObject::DATA_OBJECT()));

    void* pixvals = input->GetScalarPointer();
    int n_comp = input->GetNumberOfScalarComponents();
    if (n_comp < 1)
        vtkDebugMacro("Error reading number of components.");
    if (n_comp > 2)
        vtkDebugMacro("RGB colours detected. Using only first channel for computation.");

    std::size_t const n_points = static_cast<std::size_t>(input->GetNumberOfPoints());
    if (n_points == 0)
    {
        vtkDebugMacro("No data found!");
        return -1;
    }

    double spacing[3];
    input->GetSpacing(spacing);
    int dimensions[3];
    input->GetDimensions(dimensions);
    double origin[3];
    input->GetOrigin(origin);
    MathLib::Point3d const ll(std::array<double, 3>{{origin[0], origin[1], origin[2]}});
    GeoLib::RasterHeader const header = {
        static_cast<std::size_t>(dimensions[0]),
        static_cast<std::size_t>(dimensions[1]),
        1, ll, spacing[0], -9999};
    std::vector<double> pixels;
    pixels.reserve(n_points);
    for (std::size_t i = 0; i < n_points; ++i)
    {
        if ((n_comp == 2 || n_comp == 4) &&
            (((float*)pixvals)[(i + 1) * n_comp - 1] < 0.00000001f))
            pixels.push_back(-9999);
        else
            pixels.push_back(((float*)pixvals)[i * n_comp]);
    }
    GeoLib::Raster const* const raster(new GeoLib::Raster(header, pixels.begin(), pixels.end()));

    vtkSmartPointer<vtkPoints> new_points = vtkSmartPointer<vtkPoints>::New();
    new_points->SetNumberOfPoints(PointsPerPixel * n_points);
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    cells->Allocate(PointsPerPixel * n_points);

    double const half_cellsize(spacing[0] / 2.0);
    std::size_t pnt_idx(0);
    for (std::size_t i = 0; i < static_cast<std::size_t>(n_points); ++i)
    {
        // Skip transparent pixels
        if (n_comp == 2 || n_comp == 4)
        {
            if (((float*)pixvals)[(i + 1) * n_comp - 1] < 0.00000001f)
                continue;
        }

        double p[3];
        input->GetPoint(i, p);
        MathLib::Point3d min_pnt{std::array<double, 3>{
            {p[0] - half_cellsize, p[1] - half_cellsize, 0}}};
        MathLib::Point3d max_pnt{std::array<double, 3>{
            {p[0] + half_cellsize, p[1] + half_cellsize, 0}}};
        createPointSurface(
            new_points, cells, pnt_idx, min_pnt, max_pnt, *raster);
        pnt_idx += PointsPerPixel;
    }

    output->SetPoints(new_points);
    output->SetVerts(cells);
    output->Squeeze();

    vtkDebugMacro(<< "Created " << new_points->GetNumberOfPoints() << " points.");
    return 1;
}

void VtkImageDataToSurfacePointsFilter::createPointSurface(
    vtkSmartPointer<vtkPoints>& points,
    vtkSmartPointer<vtkCellArray>& cells,
    std::size_t pnt_idx,
    MathLib::Point3d const& min_pnt,
    MathLib::Point3d const& max_pnt,
    GeoLib::Raster const& raster)
{
    std::size_t const n_points(static_cast<std::size_t>(this->GetPointsPerPixel()));
    for (std::size_t i = 0; i < n_points; ++i)
    {
        double p[3];
        p[0] = getRandomNumber(min_pnt[0], max_pnt[0]);
        p[1] = getRandomNumber(min_pnt[1], max_pnt[1]);
        p[2] = raster.interpolateValueAtPoint(MathLib::Point3d({p[0], p[1], 0}));
        points->SetPoint(pnt_idx + i, p);
        cells->InsertNextCell(1);
        cells->InsertCellPoint(pnt_idx + i);
    }
}

double VtkImageDataToSurfacePointsFilter::getRandomNumber(double const& min, double const& max) const
{
    return ((double)std::rand() / RAND_MAX) * (max - min) + min;
}
