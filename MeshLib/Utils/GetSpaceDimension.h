/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on June 8, 2021, 12:16 PM
 */

#pragma once

#include <vector>

namespace MeshLib
{
class Node;

/**
 * \brief Computes dimension of the embedding space containing the set of given
 * points.
 *
 * The space dimension is computed by accounting the non-zero norms of
 *  \f$\mathbf x\f$, \f$\mathbf y\f$, \f$\mathbf z\f$, which are
 *  the coordinates of all nodes of the mesh. With this concept,
 *  the space dimension is:
 *    - 1, if all nodes are on a line that is parallel to
 *         \f$x\f$ axis.
 *    - 2, if all nodes are on a line that is parallel to
 *         \f$y\f$ axis.
 *    - 3, if all nodes are on a line that is parallel to
 *         \f$z\f$ axis.
 *    - 2, if all nodes are distributed on a plane that is parallel
 *         to the origin coordinate plane of \f$x-y\f$ (including the case of
 *         mesh with 1D inclined elements).
 *    - 3, if all nodes are distributed on a plane that is parallel
 *         to the origin coordinate plane of \f$y-z\f$ or \f$x-z\f$ (including
 *         the case of mesh with 1D inclined elements)..
 *    - 3, if all nodes are scattered in 3D space (e.g. mesh with inclined 1D,
 *         2D elements, 3D elements).
 */
int getSpaceDimension(std::vector<Node*> const& nodes);
};  // namespace MeshLib
