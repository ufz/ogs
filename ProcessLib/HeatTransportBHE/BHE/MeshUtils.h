// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef>
#include <unordered_map>
#include <vector>

#include "BHETypes.h"

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
    std::vector<std::vector<MeshLib::Node*>> BHE_topology_ordered_nodes;
    std::unordered_map<std::size_t, double> BHE_element_distances_from_wellhead;
    std::unordered_map<std::size_t, int> BHE_element_section_indices;

    void updateElementSectionIndices(std::vector<BHE::BHETypes> const& bhes);

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
}  // end of namespace HeatTransportBHE
}  // namespace ProcessLib
