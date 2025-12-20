// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>
#include <vector>

namespace MeshLib
{
class Element;
/**
 * \brief Element rotation matrix computation
 *
 *  This function returns a vector containing the rotation matrices of given
 * elements. The rotation matrix of an element is used to map the
 * local vector to the global coordinate system. If an element is not inclined,
 * the identity matrix is used as its rotation matrix.
 *
 * \param space_dimension The space dimension.
 * \param mesh_dimension  The mesh dimension.
 * \param elements        The mesh elements.
 * \return A vector of rotation matrices of given elements.
 */
std::vector<Eigen::MatrixXd> getElementRotationMatrices(
    int const space_dimension, int const mesh_dimension,
    std::vector<Element*> const& elements);
}  // namespace MeshLib
