/**
 * \file
 * \brief Computes mesh node ids of mesh nodes within a polygonal region, that
 * resides on the surface.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#ifdef USE_PETSC
#include <mpi.h>
#endif

#include <algorithm>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "Applications/FileIO/readGeometryFromFile.h"
#include "BaseLib/Error.h"
#include "BaseLib/FileTools.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Polygon.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshToolsLib/MeshSurfaceExtraction.h"

void writeToFile(
    std::string const& id_area_fname, std::string const& csv_fname,
    std::vector<std::pair<std::size_t, double>> const& ids_and_areas,
    std::vector<MeshLib::Node*> const& mesh_nodes)
{
    std::ofstream ids_and_area_out(id_area_fname);
    if (!ids_and_area_out)
    {
        OGS_FATAL("Unable to open the file '{:s}' - aborting.", id_area_fname);
    }
    std::ofstream csv_out(csv_fname);
    if (!csv_out)
    {
        OGS_FATAL("Unable to open the file '{:s}' - aborting.", csv_fname);
    }

    ids_and_area_out << std::setprecision(20);
    csv_out << std::setprecision(20);

    ids_and_area_out << ids_and_areas[0].first << " "
                     << ids_and_areas[0].second;
    csv_out << "ID x y z area node_id\n";  // CSV header
    csv_out << 0 << " " << *mesh_nodes[ids_and_areas[0].first]
            << ids_and_areas[0].second << " " << ids_and_areas[0].first;
    for (std::size_t k(1); k < ids_and_areas.size(); k++)
    {
        ids_and_area_out << "\n"
                         << ids_and_areas[k].first << " "
                         << ids_and_areas[k].second;
        csv_out << "\n"
                << k << " " << *mesh_nodes[ids_and_areas[k].first]
                << ids_and_areas[k].second << " " << ids_and_areas[k].first;
    }
    ids_and_area_out << "\n";
    csv_out << "\n";
}

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Computes ids of mesh nodes that are in polygonal regions and resides "
        "on the top surface. The polygonal regions have to be given in a gml- "
        "or gli-file. The found mesh nodes and the associated area are written "
        "as txt and csv data. The documentation is available at "
        "https://docs.opengeosys.org/docs/tools/model-preparation/"
        "computesurfacenodeidsinpolygonalregion.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2024, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<std::string> mesh_in(
        "m", "mesh-input-file",
        "the name of the file containing the input mesh", true, "",
        "file name of input mesh");
    cmd.add(mesh_in);
    TCLAP::ValueArg<std::string> geo_in(
        "g", "geo-file", "the name of the gml file containing the polygons",
        true, "", "file name of input geometry");
    cmd.add(geo_in);

    TCLAP::ValueArg<std::string> gmsh_path_arg("", "gmsh-path",
                                               "the path to the gmsh binary",
                                               false, "", "path as string");
    cmd.add(gmsh_path_arg);
    cmd.parse(argc, argv);

#ifdef USE_PETSC
    MPI_Init(&argc, &argv);
#endif

    std::unique_ptr<MeshLib::Mesh const> mesh(
        MeshLib::IO::readMeshFromFile(mesh_in.getValue()));
    INFO("Mesh read: {:d} nodes, {:d} elements.", mesh->getNumberOfNodes(),
         mesh->getNumberOfElements());

    GeoLib::GEOObjects geo_objs;
    FileIO::readGeometryFromFile(geo_in.getValue(), geo_objs,
                                 gmsh_path_arg.getValue());
    auto const geo_name = geo_objs.getGeometryNames()[0];
    INFO("Geometry '{:s}' read: {:d} points, {:d} polylines.",
         geo_name,
         geo_objs.getPointVec(geo_name)->size(),
         geo_objs.getPolylineVec(geo_name)->size());

    Eigen::Vector3d const dir({0.0, 0.0, -1.0});
    double angle(90);

    auto computeElementTopSurfaceAreas =
        [](MeshLib::Mesh const& mesh, Eigen::Vector3d const& d, double angle)
    {
        std::unique_ptr<MeshLib::Mesh> surface_mesh(
            MeshToolsLib::MeshSurfaceExtraction::getMeshSurface(mesh, d,
                                                                angle));
        return MeshToolsLib::MeshSurfaceExtraction::getSurfaceAreaForNodes(
            *surface_mesh);
    };

    std::vector<double> areas(computeElementTopSurfaceAreas(*mesh, dir, angle));
    std::vector<MeshLib::Node*> all_sfc_nodes(
        MeshToolsLib::MeshSurfaceExtraction::getSurfaceNodes(*mesh, dir,
                                                             angle));

    std::for_each(all_sfc_nodes.begin(), all_sfc_nodes.end(),
                  [](MeshLib::Node* p) { (*p)[2] = 0.0; });

    std::vector<MeshLib::Node*> const& mesh_nodes(mesh->getNodes());
    GeoLib::PolylineVec const* ply_vec(geo_objs.getPolylineVecObj(geo_name));
    auto const& plys(ply_vec->getVector());

    for (std::size_t j(0); j < plys.size(); j++)
    {
        if (!plys[j]->isClosed())
        {
            continue;
        }
        std::string polygon_name;
        ply_vec->getNameOfElement(plys[j], polygon_name);
        if (polygon_name.empty())
        {
            polygon_name = "Polygon-" + std::to_string(j);
        }
        // create Polygon from Polyline
        GeoLib::Polygon const polygon{*plys[j]};
        // ids of mesh nodes on surface that are within the given polygon
        std::vector<std::pair<std::size_t, double>> ids_and_areas;
        for (std::size_t k(0); k < all_sfc_nodes.size(); k++)
        {
            MeshLib::Node const& surface_node(*(all_sfc_nodes[k]));
            if (polygon.isPntInPolygon(surface_node))
            {
                ids_and_areas.emplace_back(surface_node.getID(), areas[k]);
            }
        }
        if (ids_and_areas.empty())
        {
            ERR("Polygonal part of surface '{:s}' doesn't contains nodes. No "
                "output will be generated.",
                polygon_name);
            continue;
        }

        std::string const out_path(BaseLib::extractPath(geo_in.getValue()));
        std::string id_and_area_fname(out_path + polygon_name);
        std::string csv_fname(out_path + polygon_name);
        id_and_area_fname += std::to_string(j) + ".txt";
        csv_fname += std::to_string(j) + ".csv";
        INFO(
            "Polygonal part of surface '{}' contains {} nodes. Writing to "
            "files '{}' and '{}'.",
            polygon_name,
            ids_and_areas.size(),
            id_and_area_fname,
            csv_fname);
        writeToFile(id_and_area_fname, csv_fname, ids_and_areas, mesh_nodes);
    }

    std::for_each(all_sfc_nodes.begin(), all_sfc_nodes.end(),
                  std::default_delete<MeshLib::Node>());

#ifdef USE_PETSC
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
