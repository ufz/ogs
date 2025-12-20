// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vtkType.h>
#include <vtkUnstructuredGrid.h>

#include <Eigen/Core>
#include <range/v3/view/enumerate.hpp>
#include <vector>

#include "ComputeNaturalCoordsResult.h"

namespace AU = ApplicationUtils;

/// @brief Result of an intersection of a line with a cell.
struct IntersectionResult
{
    Eigen::Vector3d point;
    double t;  // parameter along the line from p0 to p1, i.e. p = p0 + t*(p1 -
               // p0) used for sorting
};

/// @brief  Finds intersection points of a line segment with the cells of a
/// vtkUnstructuredGrid. The line segment is defined by consecutive entries in
/// realcoords.
/// @param grid The vtkUnstructuredGrid to intersect with.
/// @param realcoords original coordinates of the anchor start and end points.
/// @param free_fraction The fraction of the line segment to consider for
/// intersections.
/// @param tol Tolerance for considering points as identical.
/// @return a vector of vectors containing the intersection results which need
/// to be converted using setPhysicalPropertiesForIntersectionPoints() to fill
/// in physical properties and to obtain the final anchor format
std::vector<std::vector<IntersectionResult>> getOrderedAnchorCoords(
    vtkUnstructuredGrid* grid,
    Eigen::MatrixX3d const& realcoords,
    Eigen::VectorXd const& free_fraction,
    double const tol);

/// @brief fills the physical properties of the intersection points based on the
/// original anchor data.
/// @param anchor_coords The intersection points of the anchors with the bulk
/// mesh.
/// @param original_anchor_data The original anchor data read from the JSON
/// file.
/// @return struct containing the anchor data which will be stored in the
/// unstructured grid
AU::ComputeNaturalCoordsResult setPhysicalPropertiesForIntersectionPoints(
    std::vector<std::vector<IntersectionResult>> const& anchor_coords,
    AU::ComputeNaturalCoordsResult const& original_anchor_data);
