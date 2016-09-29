/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
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
namespace SmallDeformationWithLIE
{

void getFractureMatrixDataInMesh(
        MeshLib::Mesh const& mesh,
        std::vector<MeshLib::Element*>& vec_matrix_elements,
        std::vector<MeshLib::Element*>& vec_fracutre_elements,
        std::vector<MeshLib::Element*>& vec_fracutre_matrix_elements,
        std::vector<MeshLib::Node*>& vec_fracutre_nodes
        );

}  // namespace SmallDeformationWithLIE
}  // namespace ProcessLib
