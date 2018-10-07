/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
struct BHEMeshData
{
    std::vector<MeshLib::Element*> soil_elements;
    std::vector<int> BHE_mat_IDs;
    std::vector<std::vector<MeshLib::Element*>> BHE_elements;
    std::vector<MeshLib::Node*> soil_nodes;
    std::vector<std::vector<MeshLib::Node*>> BHE_nodes;
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
