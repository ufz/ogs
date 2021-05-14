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

std::vector<Eigen::MatrixXd> getElementRotationMatrices(
    int const space_dimension, int const mesh_dimension,
    std::vector<Element*> const& elements);
}  // namespace MeshLib
