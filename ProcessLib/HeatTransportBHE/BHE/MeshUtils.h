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
/**
 * get data about fracture and matrix elements/nodes from a mesh
 *
 * @param mesh  A mesh which includes BHE elements, i.e. 1-dimensional
 * elements. It is assumed that elements forming a BHE have a distinct
 * material ID.
 * @param vec_soil_elements a vector of soil elements
 * @param vec_BHE_elements  a vector of BHE elements (grouped by BHE IDs)
 * @param vec_BHE_nodes   a vector of BHE nodes (grouped by BHE IDs)
 */
void getBHEDataInMesh(
    MeshLib::Mesh const& mesh,
    std::vector<MeshLib::Element*>& vec_soil_elements,
    std::vector<int>& vec_BHE_mat_IDs,
    std::vector<std::vector<MeshLib::Element*>>& vec_BHE_elements,
    std::vector<MeshLib::Node*>& vec_pure_soil_nodes,
    std::vector<std::vector<MeshLib::Node*>>& vec_BHE_nodes);
}  // end of namespace HeatTransportBHE
}  // namespace ProcessLib
