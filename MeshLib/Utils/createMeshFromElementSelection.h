// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
