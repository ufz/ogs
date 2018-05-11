/**
 * \file
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vtkPolyDataAlgorithm.h>

#include "GeoLib/Raster.h"
#include "VtkAlgorithmProperties.h"

class vtkInformation;
class vtkInformationVector;

/// A VTK filter that creates a point cloud representing a surface defined by
/// pixel values
class VtkImageDataToSurfacePointsFilter : public vtkPolyDataAlgorithm
{
public:
    /// Create a new objects (required because of VTKs reference counting)
    static VtkImageDataToSurfacePointsFilter* New();

    vtkTypeMacro(VtkImageDataToSurfacePointsFilter, vtkPolyDataAlgorithm);

    /// @brief Prints information about itself.
    void PrintSelf(ostream& os, vtkIndent indent) override;

    vtkGetMacro(PointsPerPixel, vtkIdType);
    vtkSetMacro(PointsPerPixel, vtkIdType);

protected:
    VtkImageDataToSurfacePointsFilter();
    ~VtkImageDataToSurfacePointsFilter() override;

    /// Sets input port to vtkImageData.
    int FillInputPortInformation(int port, vtkInformation* info) override;

    /// Updates the graphical object
    int RequestData(vtkInformation* request, vtkInformationVector** inputVector,
                    vtkInformationVector* outputVector) override;

private:
    VtkImageDataToSurfacePointsFilter(const VtkImageDataToSurfacePointsFilter&) = delete;
    void operator=(const VtkImageDataToSurfacePointsFilter&) = delete;

    /// Returns a random number in [min, max]
    void createPointSurface(vtkSmartPointer<vtkPoints>& points,
                            vtkSmartPointer<vtkCellArray>& cells,
                            std::size_t pnt_idx,
                            MathLib::Point3d const& min_pnt,
                            MathLib::Point3d const& max_pnt,
                            GeoLib::Raster const& raster);

    /// Returns a random number in [min, max]
    double getRandomNumber(double const& min, double const& max) const;

    vtkIdType PointsPerPixel;
};
