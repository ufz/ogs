/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>
#include <memory>

namespace GeoLib
{
class Polyline;
class Surface;
class GEOObjects;
}  // namespace GeoLib

namespace FileIO
{
/// Creates a plane surface from the given polyline. The polyline has to be
/// closed, i.e. the first and the last point have to be the identical. The
/// triangulation of the polyline is done by the finite element meshing tool
/// Gmsh. Finally, the resulting mesh is converted into a GeoLib::Surface which
/// is inserted into the \c GeoLib::GEOObjects instance \c geometries using the
/// name \c geometry_name.
bool createSurface(GeoLib::Polyline const& ply,
                   GeoLib::GEOObjects& geometries,
                   std::string const& geometry_name,
                   std::string const& gmsh_binary);

/// Creates a plane surface from the given polyline. The polyline has to be
/// closed, i.e. the first and the last point have to be the identical. The
/// triangulation of the polyline is done via a simple ear clipping algorithm.
std::unique_ptr<GeoLib::Surface> createSurfaceWithEarClipping(
    GeoLib::Polyline const& line);
}  // namespace FileIO
