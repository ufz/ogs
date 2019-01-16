
/**
 *
 * @copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#pragma once

#include <string>

#include "MeshLib/Mesh.h"

#include "GocadSGridReader.h"

namespace FileIO
{
namespace Gocad
{

void generateFaceSets(GocadSGridReader const& reader, std::string const& path);

}  //  namespace Gocad
}  //  namespace FileIO
