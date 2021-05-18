/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on May 14, 2021, 2:38 PM
 */

#pragma once

#include <Eigen/Dense>
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
 * @param space_dimension The space dimension.
 * @param mesh_dimension  The mesh dimension.
 * @param elements        The mesh elements.
 * @return A vector of rotation matrices of given elements.
 */
std::vector<Eigen::MatrixXd> getElementRotationMatrices(
    int const space_dimension, int const mesh_dimension,
    std::vector<Element*> const& elements);
}  // namespace MeshLib
