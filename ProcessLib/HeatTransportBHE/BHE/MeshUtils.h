/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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
 * @param mesh  A mesh which includes BHE elements, i.e. 1-dimensional
 * elements. It is assumed that elements forming a BHE have a distinct
 * material ID.
 */
BHEMeshData getBHEDataInMesh(MeshLib::Mesh const& mesh);
}  // end of namespace HeatTransportBHE
}  // namespace ProcessLib
