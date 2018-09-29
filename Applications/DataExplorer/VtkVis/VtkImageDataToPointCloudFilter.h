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

#include "MathLib/Point3d.h"
#include "VtkAlgorithmProperties.h"

class vtkInformation;
class vtkInformationVector;

/// A VTK Filter that will create a point cloud with local densities based on
/// pixel values.
class VtkImageDataToPointCloudFilter : public vtkPolyDataAlgorithm
{
public:
    /// Create a new objects (required because of VTKs reference counting)
    static VtkImageDataToPointCloudFilter* New();

    vtkTypeMacro(VtkImageDataToPointCloudFilter, vtkPolyDataAlgorithm);

    /// Prints information about the filter
    void PrintSelf(ostream& os, vtkIndent indent) override;

    void useLinearInterpolation() { SetIsLinear(true); }
    void useLogarithmicInterpolation(double gamma)
    {
        SetIsLinear(false);
        SetGamma(gamma);
    }

    vtkGetMacro(IsLinear, bool);
    vtkSetMacro(IsLinear, bool);

    vtkGetMacro(MinNumberOfPointsPerCell, vtkIdType);
    vtkSetMacro(MinNumberOfPointsPerCell, vtkIdType);

    vtkGetMacro(MaxNumberOfPointsPerCell, vtkIdType);
    vtkSetMacro(MaxNumberOfPointsPerCell, vtkIdType);

    vtkGetMacro(Gamma, double);
    vtkSetMacro(Gamma, double);

    vtkGetMacro(PointScaleFactor, double);
    vtkSetMacro(PointScaleFactor, double);

    vtkGetMacro(MinHeight, double);
    vtkSetMacro(MinHeight, double);

    vtkGetMacro(MaxHeight, double);
    vtkSetMacro(MaxHeight, double);

protected:
    VtkImageDataToPointCloudFilter();
    ~VtkImageDataToPointCloudFilter() = default;

    /// Sets input port to vtkImageData.
    int FillInputPortInformation(int port, vtkInformation* info) override;

    /// Updates the graphical object
    int RequestData(vtkInformation* request, vtkInformationVector** inputVector,
                    vtkInformationVector* outputVector) override;

    double Gamma;
    double PointScaleFactor;
    double MinHeight;
    double MaxHeight;
    vtkIdType MinNumberOfPointsPerCell;
    vtkIdType MaxNumberOfPointsPerCell;
    bool IsLinear;


private:
    /// Creates the point objects based on the parameters set by the user
    void createPoints(vtkSmartPointer<vtkPoints>& points,
                      vtkSmartPointer<vtkCellArray>& cells, std::size_t pnt_idx,
                      std::size_t n_points, MathLib::Point3d const& min_pnt,
                      MathLib::Point3d const& max_pnt) const;

    /// Returns a random number in [min, max]
    double getRandomNumber(double const& min, double const& max) const;

    /// Interpolates the required number of points. Using gamma = 1 this is a linear
    /// interpolation, for values smaller/greater than 1 the curve becomes
    /// logarithmic/exponential, resulting in a smaller/larger difference of point
    /// densities given the same input values.
    std::size_t interpolate(double low, double high, double p, double gamma) const;
};
