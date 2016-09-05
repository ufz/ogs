/*
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef CONVERTMESHNODESTOGEOMETRY_H
#define CONVERTMESHNODESTOGEOMETRY_H

#include <memory>
#include <string>
#include <vector>

#include "GeoLib/GEOObjects.h"
#include "GeoLib/Point.h"

namespace MeshLib
{
class Mesh;
}

namespace MeshGeoToolsLib
{

void convertMeshNodesToGeometry(MeshLib::Mesh const& mesh,
                                std::vector<std::size_t> const& node_ids,
                                std::string& geo_name,
                                GeoLib::GEOObjects& geometry_sets);

}  // end namespace MeshGeoToolsLib

#endif
