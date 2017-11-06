/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <Eigen/Eigen>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"

namespace ProcessLib
{
namespace LIE
{
/// Compute a normal vector of the given element
///
/// Computed normal vector is oriented in the left direction of the given line
/// element such that computeRotationMatrix() returns the indentity matrix for
/// line elements parallel to a vector (1,0,0)
void computeNormalVector(MeshLib::Element const& e, unsigned const global_dim,
                         Eigen::Vector3d& normal_vector);

/// Compute a rotation matrix from global to local coordinates using the given
/// elements' normal vector and based on the two vectors forming the element's
/// surface plane.
///
/// In the 2D case (line element) the resulting y axis should be same as the
/// given normal vector. In the 3D case (tri, quad, etc.) the resulting z-axis
/// should be the same as the given normal vector.
///
/// \param e the element.
/// \param n the element's normal.
/// \param global_dim the space dimension in which the element is embedded.
/// \param R the output rotation matrix.
void computeRotationMatrix(MeshLib::Element const& e, Eigen::Vector3d const& n,
                           unsigned const global_dim, Eigen::MatrixXd& R);

/// compute physical coordinates from the given shape vector, i.e. from the
/// natural coordinates
template <typename Derived>
MathLib::Point3d computePhysicalCoordinates(
    MeshLib::Element const& e, Eigen::MatrixBase<Derived> const& shape)
{
    MathLib::Point3d pt;
    for (unsigned i = 0; i < e.getNumberOfNodes(); i++)
    {
        MeshLib::Node const& node = *e.getNode(i);
        for (unsigned j = 0; j < 3; j++)
            pt[j] += shape[i] * node[j];
    }
    return pt;
}

}  // namespace LIE
}  // namespace ProcessLib
