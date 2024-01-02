/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

namespace MeshLib
{
class Mesh;
class Element;

/// Creates a new mesh from a vector of elements.
///
/// \note The elements are owned by the returned mesh object as well as the
/// nodes and will be destructed together with the mesh.
std::unique_ptr<Mesh> createMeshFromElementSelection(
    std::string mesh_name, std::vector<Element*> const& elements);
}  // namespace MeshLib
