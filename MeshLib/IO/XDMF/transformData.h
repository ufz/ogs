/**
 * \file
 * \author Tobias Meisel
 * \date   2020-11-13
 * \brief  Transforms OGS Mesh into vectorized data
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <vector>

#include "XdmfData.h"

namespace MeshLib
{
class Mesh;
}

namespace MeshLib::IO
{
std::vector<AttributeMeta> transformAttributes(MeshLib::Mesh const& mesh);
Geometry transformGeometry(MeshLib::Mesh const& mesh);
Topology transformTopology(MeshLib::Mesh const& mesh);
}  // namespace MeshLib::IO