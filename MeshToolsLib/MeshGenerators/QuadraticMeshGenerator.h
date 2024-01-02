/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <memory>

namespace MeshLib
{
class Mesh;
}

namespace MeshToolsLib
{
/// Create a quadratic order mesh from the linear order mesh. For some element
/// types like Quad-4, a centre node might be added if the \c add_centre_node
/// flag is set, yielding a Quad-9.
std::unique_ptr<MeshLib::Mesh> createQuadraticOrderMesh(
    MeshLib::Mesh const& linear_mesh, bool const add_centre_node);

}  // namespace MeshToolsLib
