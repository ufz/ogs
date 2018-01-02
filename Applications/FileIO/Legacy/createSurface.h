/**
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace GeoLib
{
class Polyline;
class GEOObjects;
}

namespace FileIO
{
/// Creates a plane surface from the given polyline. The polyline have to be
/// closed, i.e. the first and the last point have to be the identical. The
/// triangulation of the polyline is done by the finite element meshing tool
/// Gmsh. Finally, the resulting mesh is converted into a GeoLib::Surface which
/// is inserted into the \c GeoLib::GEOObjects instance \c geometries using the
/// name \c geometry_name.
bool createSurface(GeoLib::Polyline const& polyline,
                   GeoLib::GEOObjects& geometries,
                   std::string const& geometry_name);
}
