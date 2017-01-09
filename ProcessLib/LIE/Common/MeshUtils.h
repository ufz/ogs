/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <vector>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

namespace ProcessLib
{
namespace LIE
{

/**
 * get data about fracture and matrix elements/nodes from a mesh
 *
 * @param mesh  A mesh which includes fracture elements, i.e. lower-dimensional elements.
 * It is assumed that elements forming a fracture have a distinct material ID.
 * @param vec_matrix_elements  a vector of matrix elements
 * @param vec_fracture_mat_IDs  fracture material IDs found in the mesh
 * @param vec_fracture_elements  a vector of fracture elements (grouped by fracture IDs)
 * @param vec_fracture_matrix_elements  a vector of fracture elements and matrix elements
 * connecting to the fracture (grouped by fracture IDs)
 * @param vec_fracture_nodes  a vector of fracture element nodes (grouped by fracture IDs)
 */
void getFractureMatrixDataInMesh(
        MeshLib::Mesh const& mesh,
        std::vector<MeshLib::Element*>& vec_matrix_elements,
        std::vector<int>& vec_fracture_mat_IDs,
        std::vector<std::vector<MeshLib::Element*>>& vec_fracture_elements,
        std::vector<std::vector<MeshLib::Element*>>& vec_fracture_matrix_elements,
        std::vector<std::vector<MeshLib::Node*>>& vec_fracture_nodes
        );

}  // namespace LIE
}  // namespace ProcessLib
