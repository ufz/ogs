/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHGENERATOR_H_
#define MESHGENERATOR_H_

#include "GeoLib/Point.h"
#include "MeshLib/Mesh.h"

namespace MeshLib
{
namespace MeshGenerator
{
/**
 * Generate an 1D Line-Element mesh. The mesh is generated in x-direction.
 *
 * \param length Mesh's length in x-direction.
 * \param subdivision Number of subdivisions.
 * \param origin Optional mesh's origin (the most left point) with GeoLib::ORIGIN default.
 */
Mesh* generateLineMesh(const double length,
                       const std::size_t subdivision,
                       GeoLib::Point const& origin = GeoLib::ORIGIN);

/**
 * Generate a regular 2D Quad-Element mesh. The mesh is generated in the
 * x-y-plane.
 *
 * \param length Mesh's dimensions in x- and y-directions.
 * \param subdivision Number of subdivisions.
 * \param origin Optional mesh's origin (the lower left corner) with GeoLib::ORIGIN default.
 */
Mesh* generateRegularQuadMesh(const double length,
                              const std::size_t subdivision,
                              GeoLib::Point const& origin = GeoLib::ORIGIN);

/**
 * Generate a regular 3D Hex-Element mesh.
 *
 * \param length Mesh's dimensions in x- and y- and z-directions.
 * \param subdivision Number of subdivisions.
 * \param origin Optional mesh's origin (the bottom lower left corner) with GeoLib::ORIGIN default.
 */
Mesh* generateRegularHexMesh(const double length,
                              const std::size_t subdivision,
                              GeoLib::Point const& origin = GeoLib::ORIGIN);
}  //MeshGenerator
} //MeshLib

#endif //MESHGENERATOR_H_
