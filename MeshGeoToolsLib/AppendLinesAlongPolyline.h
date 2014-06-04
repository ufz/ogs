/**
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifndef APPENDLINESALONGPOLYLINES_H_
#define APPENDLINESALONGPOLYLINES_H_

// GeoLib
#include "PolylineVec.h"

namespace MeshLib
{
class Mesh;
}

namespace MeshGeoToolsLib
{
MeshLib::Mesh* appendLinesAlongPolylines(const MeshLib::Mesh &mesh, const GeoLib::PolylineVec &ply_vec);
}

#endif
