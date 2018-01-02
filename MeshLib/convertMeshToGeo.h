/**
 * \file
 * \author Karsten Rink
 * \date   2013-07-05
 * \brief  Definition of mesh to geometry conversion.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <limits>
#include <string>

namespace GeoLib
{
class GEOObjects;
class Surface;
}

namespace MeshLib
{

    class Mesh;

    /**
     * Converts a 2D mesh into a geometry.
     * A new geometry with the name of the mesh will be inserted into geo_objects, consisting
     * of points identical with mesh nodes and one surface representing the mesh. Triangles are
     * converted to geometric triangles, quads are split into two triangles, all other elements
     * are ignored.
     */
    bool convertMeshToGeo(const MeshLib::Mesh &mesh, GeoLib::GEOObjects &geo_objects, double eps = std::numeric_limits<double>::epsilon());

    /**
     * Converts a surface into a triangular mesh
     * @param sfc         Surface object
     * @param mesh_name   New mesh name
     * @param eps         Minimum distance for nodes not to be collapsed
     * @return a pointer to a converted mesh object. nullptr is returned if the conversion fails.
     */
    MeshLib::Mesh* convertSurfaceToMesh(const GeoLib::Surface &sfc, const std::string &mesh_name, double eps = std::numeric_limits<double>::epsilon());

} // namespace
