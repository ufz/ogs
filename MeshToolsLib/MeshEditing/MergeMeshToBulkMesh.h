// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <string>
#include <unordered_map>

namespace MeshLib
{
class Mesh;
class Node;
class Element;
}  // namespace MeshLib

namespace MeshToolsLib
{
/**
 * Merges the mesh \c other_mesh into the bulk mesh \c bulk_mesh.
 * Optionally sets constant integration point stress and primary nodal variables
 * in the bulk mesh to specified values for the merged part.
 *
 * Note: This function moves nodes and elements from \c other_mesh to the merged
 *       mesh (\c bulk_mesh). It deletes the duplicate nodes of \c other_mesh at
 *       the mesh interface. The merged mesh manages memory deallocation for
 *       nodes and elements. The member function `shallowClean` must be called
 *       on \c bulk_mesh and \c other_mesh afterward to prevent memory leaks.
 *
 * @param bulk_mesh           Bulk mesh.
 * @param other_mesh          Mesh to be merged.
 * @param initial_value_dict  Initial values to be specified for the merged
 *                            part.
 * @return
 */
std::unique_ptr<MeshLib::Mesh> mergeMeshToBulkMesh(
    MeshLib::Mesh const& bulk_mesh, MeshLib::Mesh const& other_mesh,
    std::unordered_map<std::string, double>& initial_value_dict);

}  // namespace MeshToolsLib
