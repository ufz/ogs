/**
 * \file
 * \author Karsten Rink
 * \date   2012-09-27
 * \brief  Definition of readMeshFromFile function.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <string>

namespace MeshLib
{
class Mesh;

namespace IO
{
MeshLib::Mesh* readMeshFromFile(const std::string& file_name,
                                bool const compute_element_neighbors = false);
}
}
