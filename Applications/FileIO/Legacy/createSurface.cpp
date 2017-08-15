/**
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cstdio>
#include <list>
#include <memory>

#include <logog/include/logog.hpp>

#include "createSurface.h"

#include "Applications/FileIO/Gmsh/GMSHInterface.h"
#include "Applications/FileIO/Gmsh/GmshReader.h"

#include "GeoLib/GEOObjects.h"
#include "GeoLib/Point.h"
#include "GeoLib/Polygon.h"
#include "GeoLib/Polyline.h"
#include "GeoLib/Surface.h"

#include "MeshLib/convertMeshToGeo.h"
#include "MeshLib/Mesh.h"

namespace FileIO
{
bool createSurface(GeoLib::Polyline const& ply,
                   GeoLib::GEOObjects& geometries,
                   std::string const& geometry_name)
{
    if (!ply.isClosed())
    {
        WARN("Error in createSurface() - Polyline is not closed.");
        return false;
    }

    if (ply.getNumberOfPoints() <= 2)
    {
        WARN(
            "Error in createSurface() - Polyline consists of less "
            "than three points and therefore cannot be triangulated.");
        return false;
    }

    // create new GEOObjects and insert a copy of the polyline
    auto polyline_points = std::make_unique<std::vector<GeoLib::Point*>>();
    GeoLib::GEOObjects geo;
    auto ply_points = ply.getPointsVec();
    for (auto p : ply_points)
        polyline_points->push_back(new GeoLib::Point(*p));
    std::string ply_name = "temporary_polyline_name";
    geo.addPointVec(std::move(polyline_points), ply_name);
    auto polyline =
        std::make_unique<GeoLib::Polyline>(*geo.getPointVec(ply_name));
    for (std::size_t k(0); k < ply.getNumberOfPoints(); ++k)
        polyline->addPoint(ply.getPointID(k));
    auto polylines = std::make_unique<std::vector<GeoLib::Polyline*>>();
    polylines->push_back(polyline.release());
    geo.addPolylineVec(std::move(polylines), ply_name);

    // use GMSHInterface to create a mesh from the closed polyline
    std::vector<std::string> geo_names;
    geo.getGeometryNames(geo_names);
    FileIO::GMSH::GMSHInterface gmsh_io(
        geo, false, FileIO::GMSH::MeshDensityAlgorithm::FixedMeshDensity, 0.0,
        0.0, 0.0, geo_names, false, false);
    gmsh_io.setPrecision(std::numeric_limits<double>::digits10);

    char file_base_name_c[L_tmpnam];
    if (! std::tmpnam(file_base_name_c))
    {
       OGS_FATAL("Could not create unique file name.");
    }
    std::string const file_base_name(file_base_name_c);
    gmsh_io.writeToFile(file_base_name + ".geo");
    std::string gmsh_command =
        "gmsh -2 -algo meshadapt \"" + file_base_name + ".geo\"";
    gmsh_command += " -o \"" + file_base_name + ".msh\"";
    std::system(gmsh_command.c_str());
    auto surface_mesh =
        FileIO::GMSH::readGMSHMesh("\"" + file_base_name + ".msh\"");
    if (!surface_mesh)
    {
        WARN("The surface mesh could not be created.");
        return false;
    }
    if (std::remove((file_base_name + ".geo").c_str()) !=0)
        WARN("Could not remove temporary file '%s'.",
            (file_base_name + ".geo").c_str());
    if (std::remove((file_base_name + ".msh").c_str()) !=0)
        WARN("Could not remove temporary file '%s'.",
            (file_base_name + ".msh").c_str());

    // convert the surface mesh into a geometric surface
    if (!MeshLib::convertMeshToGeo(*surface_mesh, geometries,
                                   std::numeric_limits<double>::epsilon()))
    {
        WARN("The surface mesh could not be converted to a geometry.");
        return false;
    }
    std::string merged_geometry_name("geometry_with_surfaces");
    geometries.mergeGeometries({geometry_name, surface_mesh->getName()},
                               merged_geometry_name);
    geometries.removeSurfaceVec(geometry_name);
    geometries.removePolylineVec(geometry_name);
    geometries.removePointVec(geometry_name);
    geometries.removeSurfaceVec(surface_mesh->getName());
    geometries.removePolylineVec(surface_mesh->getName());
    geometries.removePointVec(surface_mesh->getName());
    geometries.renameGeometry(merged_geometry_name, geometry_name);

    return true;
}

} // end namespace
