// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

namespace MeshLib
{
class Element;
class Mesh;
class Node;
}  // namespace MeshLib

namespace ProcessLib
{
namespace HeatTransportBHE
{
/* TODO (naumov) Just an idea
struct BheMeshSubset
{
    int material_id;
    std::vector<MeshLib::Element*> elements;
    std::vector<MeshLib::Node*> nodes;
};
*/

struct BHEMeshData
{
    std::vector<int> BHE_mat_IDs;
    std::vector<std::vector<MeshLib::Element*>> BHE_elements;
    std::vector<std::vector<MeshLib::Node*>> BHE_nodes;

    // TODO (naumov) Just an idea: std::vector<BheMeshSubset> mesh_subsets;
};

/**
 * get data about fracture and matrix elements/nodes from a mesh
 *
 * \param mesh  A mesh which includes BHE elements, i.e. 1-dimensional
 * elements. It is assumed that elements forming a BHE have a distinct
 * material ID.
 */
BHEMeshData getBHEDataInMesh(MeshLib::Mesh const& mesh);

/// Inlet and outlet endpoint nodes of a BHE line-element chain, derived
/// from the elements' local node ordering. The inlet is the chain endpoint
/// that appears as node 0 of its (only) connected line element; the outlet
/// is the other endpoint. This is the same node-0 -> node-1 convention used
/// by the local assembler when building `_element_direction = (p1 - p0)`.
struct BHEEndpoints
{
    MeshLib::Node const* inlet;
    MeshLib::Node const* outlet;
};

/// Derive the inlet/outlet endpoint nodes of a BHE from the line-element
/// chain's local node ordering. Also validates that the elements form a
/// single connected chain with a consistent orientation (every element's
/// node 0 equals the previous element's node 1); calls OGS_FATAL with a
/// diagnostic message naming the offending node on inconsistency.
BHEEndpoints findBHEEndpointsFromElementOrdering(
    std::vector<MeshLib::Element*> const& bhe_elements);

}  // end of namespace HeatTransportBHE
}  // namespace ProcessLib
