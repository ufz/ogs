/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>
#include <vector>

namespace MeshLib
{
class Element;
}

namespace ApplicationUtils
{
/// Write elements as METIS graph file
/// \param elements The mesh elements.
/// \param file_name File name with an extension of mesh.
void writeMETIS(std::vector<MeshLib::Element*> const& elements,
                const std::string& file_name);

}  // namespace ApplicationUtils
