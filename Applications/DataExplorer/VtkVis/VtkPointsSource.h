// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

// ** INCLUDES **
#include "VtkAlgorithmProperties.h"
#include <vtkPolyDataAlgorithm.h>

#include "GeoLib/Point.h"

/**
 * \brief VtkPointsSource is a VTK source object for the visualization
 * of point data. As a vtkPolyDataAlgorithm it outputs polygonal data.
 */
class VtkPointsSource : public vtkPolyDataAlgorithm, public VtkAlgorithmProperties
{
public:
    /// Create new objects with New() because of VTKs object reference counting.
    static VtkPointsSource* New();

    vtkTypeMacro(VtkPointsSource,vtkPolyDataAlgorithm);

    /// Sets the points as a vector
    void setPoints(const std::vector<GeoLib::Point*>* points) { _points = points; }

    /// Prints its data on a stream.
    void PrintSelf(ostream& os, vtkIndent indent) override;

    void SetUserProperty(QString name, QVariant value) override;

protected:
    VtkPointsSource();
    ~VtkPointsSource() override = default;

    /// Computes the polygonal data object.
    int RequestData(vtkInformation* request,
                    vtkInformationVector** inputVector,
                    vtkInformationVector* outputVector) override;

    int RequestInformation(vtkInformation* request,
                           vtkInformationVector** inputVector,
                           vtkInformationVector* outputVector) override;

    /// The points to visualize
    const std::vector<GeoLib::Point*>* _points{nullptr};

private:
};
