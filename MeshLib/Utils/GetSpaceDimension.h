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
 *  the space dimension of a mesh is:
 *    - 1, if the mesh is 1D and and the mesh is parallel either to
 *         \f$x\f$, \f$y\f$ or to \f$z\f$ axis.
 *    - 2, if the mesh is 1D but it contains inclined elements on the origin
 *         coordinate plane of \f$x-y,\, y-z,\, \text{or }\, x-z\f$.
 *         That means the coordinates of all nodes are
 *         \f$(x, y, 0)\f$, \f$(x, 0, z)\f$ or \f$(0, y, z)\f$.
 *    - 3, if the mesh is 1D and but it is not on any origin
 *         coordinate plane of \f$x-y,\, y-z,\, \text{or }\, x-z\f$.
 *    - 2, if the mesh is 2D and it is on the origin
 *         coordinate plane of \f$x-y,\, y-z,\, \text{or }\, x-z\f$.
 *    - 2, if the mesh is 2D and it is on vertical or horizontal plane
 *         that is parallel to the the origin coordinate
 *         plane of \f$x-y,\, y-z,\, \text{or }\, x-z\f$ but with an offset.
 *    - 3, if the mesh contains inclined 2D elements.
 *    - 3, if the mesh contains 3D elements.
 */
int getSpaceDimension(std::vector<Node*> const& nodes);
};  // namespace MeshLib
