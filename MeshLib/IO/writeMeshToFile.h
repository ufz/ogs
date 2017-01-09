/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
int writeMeshToFile(const MeshLib::Mesh &mesh, const std::string &file_name);
}
}
