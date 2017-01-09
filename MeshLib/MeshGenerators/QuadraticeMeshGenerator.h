/**
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 */

#include <memory>

#include "MeshLib/Mesh.h"

namespace MeshLib
{

/// create a quadratic order mesh from the linear order mesh
std::unique_ptr<Mesh> createQuadraticOrderMesh(Mesh const& linear_order_mesh);

} // namespace MeshLib
