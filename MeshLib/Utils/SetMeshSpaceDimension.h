/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on May 5, 2021, 12:46 PM
 */

#pragma once

#include <memory>
#include <vector>

namespace MeshLib
{
class Mesh;

/**
 * In this function, the dimension of the space, which contains all nodes of the
 * bulk mesh, is computed, and it is then set as a member of the elements of all
 * meshes.
 *
 * The space dimension is needed for the numerical simulations with a bulk mesh
 * with elements with different dimensions. For example, in a bulk mesh for the
 * numerical simulation of the liquid flow in the fractured porous media, the 1D
 * or 2D elements can be used to represent the fractures, while the elements
 * with higher dimension can be used to discretise the porous medium matrix
 * domain.
 *
 * @param meshes All meshes include the bulk mesh and the meshes for boundary
 * conditions.
 */
void setMeshSpaceDimension(std::vector<std::unique_ptr<Mesh>> const& meshes);
};  // namespace MeshLib
