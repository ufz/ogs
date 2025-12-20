// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vtkType.h>

#include <Eigen/Core>

namespace ApplicationUtils
{
struct ComputeNaturalCoordsResult
{
    Eigen::MatrixXd natural_coords;
    Eigen::MatrixXd real_coords;
    Eigen::VectorXd initial_anchor_stress;
    Eigen::VectorXd maximum_anchor_stress;
    Eigen::VectorXd residual_anchor_stress;
    Eigen::VectorXd anchor_cross_sectional_area;
    Eigen::VectorXd anchor_stiffness;
    Eigen::VectorX<vtkIdType> bulk_element_ids;
    Eigen::VectorX<vtkIdType> point_cloud_node_ids;
    bool success;
};
}  // namespace ApplicationUtils
