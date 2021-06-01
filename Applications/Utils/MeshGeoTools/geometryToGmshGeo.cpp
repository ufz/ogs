/**
 * \file
 * \brief  A small tool to create a gmsh geometry out of gml geometry.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

// ThirdParty
#include <tclap/CmdLine.h>

#include "Applications/FileIO/Gmsh/GMSHInterface.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/IO/XmlIO/Boost/BoostXmlGmlInterface.h"
#include "InfoLib/GitInfo.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Tool to create gmsh geometry (geo-file) out of a gml geometry.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2021, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<std::string> geo_output_arg(
        "o", "output", "output gmsh geometry file (*.geo)", true, "", "output file");
    cmd.add(geo_output_arg);
    TCLAP::MultiArg<std::string> geo_input_arg(
        "i", "input", "input geometry file (*.gml)", true, "input file name");
    cmd.add(geo_input_arg);
    TCLAP::ValueArg<unsigned> max_number_of_points_in_quadtree_leaf_arg(
        "", "max_points_in_quadtree_leaf", "positive number", false, 2,
        "max points in a quadtree leaf before the leaf is split");
    cmd.add(max_number_of_points_in_quadtree_leaf_arg);
    TCLAP::ValueArg<double> mesh_density_scaling_points_arg(
        "", "mesh_density_scaling_at_points", "positive floating point number",
        false, 0.2, "desired mesh density at points");
    cmd.add(mesh_density_scaling_points_arg);
    TCLAP::ValueArg<double> mesh_density_scaling_stations_arg(
        "", "mesh_density_scaling_at_stations",
        "positive floating point number", false, 0.05,
        "desired mesh density at stations");
    cmd.add(mesh_density_scaling_stations_arg);
    TCLAP::ValueArg<double> average_point_density_arg(
        "a", "average_point_density",
        "average point density / average edge length as a positive floating "
        "point number",
        false, 1,
        "desired point density / edge length in homogeneous meshing approach");
    cmd.add(average_point_density_arg);
    TCLAP::SwitchArg homogeneous_flag(
        "", "homogeneous", "Use Gmsh homogeneous meshing method.", false);
    cmd.add(homogeneous_flag);

    cmd.parse(argc, argv);

    GeoLib::GEOObjects geo_objects;
    for (auto const& geometry_name : geo_input_arg.getValue())
    {
        GeoLib::IO::BoostXmlGmlInterface xml(geo_objects);
        try
        {
            if (!xml.readFile(geometry_name))
            {
                return EXIT_FAILURE;
            }
        }
        catch (std::runtime_error const& err)
        {
            ERR("Failed to read file '{:s}'.", geometry_name);
            ERR("{:s}", err.what());
            return EXIT_FAILURE;
        }
        INFO("Successfully read file '{:s}'.", geometry_name);
    }

    auto const geo_names = geo_objects.getGeometryNames();

    bool const rotate = false;
    bool const keep_preprocessed_geometry = false;

    if (homogeneous_flag.getValue())
    {  // homogeneous meshing
        double const average_mesh_density =
            average_point_density_arg.getValue();
        FileIO::GMSH::GMSHInterface gmsh_io(
            geo_objects, true,
            FileIO::GMSH::MeshDensityAlgorithm::FixedMeshDensity,
            average_mesh_density, 0.0, 0, geo_names, rotate,
            keep_preprocessed_geometry);
        BaseLib::IO::writeStringToFile(gmsh_io.writeToString(),
                                       geo_output_arg.getValue());
    }
    else
    {  // adaptive meshing
        unsigned const max_number_of_points_in_quadtree_leaf =
            max_number_of_points_in_quadtree_leaf_arg.getValue();
        double const mesh_density_scaling_points =
            mesh_density_scaling_points_arg.getValue();
        double const mesh_density_scaling_stations =
            mesh_density_scaling_stations_arg.getValue();
        FileIO::GMSH::GMSHInterface gmsh_io(
            geo_objects, true,
            FileIO::GMSH::MeshDensityAlgorithm::AdaptiveMeshDensity,
            mesh_density_scaling_points, mesh_density_scaling_stations,
            max_number_of_points_in_quadtree_leaf, geo_names, rotate,
            keep_preprocessed_geometry);
        BaseLib::IO::writeStringToFile(gmsh_io.writeToString(),
                                       geo_output_arg.getValue());
    }

    return EXIT_SUCCESS;
}
