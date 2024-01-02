/**
 * \file
 * 2012/03/07 KR Initial implementation
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include <fstream>
#include <memory>
#include <string>

#ifdef USE_PETSC
#include <mpi.h>
#endif

#include "BaseLib/FileTools.h"
#include "GeoLib/AABB.h"
#include "GeoLib/IO/AsciiRasterInterface.h"
#include "GeoLib/Raster.h"
#include "InfoLib/GitInfo.h"
#include "MathLib/MathTools.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSearch/MeshElementGrid.h"
#include "MeshLib/Node.h"
#include "MeshToolsLib/MeshEditing/ProjectPointOnMesh.h"
#include "MeshToolsLib/MeshGenerators/MeshLayerMapper.h"

double getClosestPointElevation(MeshLib::Node const& p,
                                std::vector<MeshLib::Node*> const& nodes,
                                double const& max_dist)
{
    double sqr_shortest_dist(max_dist);
    double elevation(p[2]);
    for (MeshLib::Node* node : nodes)
    {
        double sqr_dist = (p[0] - (*node)[0]) * (p[0] - (*node)[0]) +
                          (p[1] - (*node)[1]) * (p[1] - (*node)[1]);
        if (sqr_dist < sqr_shortest_dist)
        {
            sqr_shortest_dist = sqr_dist;
            elevation = (*node)[2];
        }
    }
    return elevation;
}

int main(int argc, char* argv[])
{
#ifdef USE_PETSC
    MPI_Init(&argc, &argv);
#endif
    TCLAP::CmdLine cmd(
        "Changes the elevation of 2D mesh nodes based on either raster data or "
        "another 2D mesh. In addition, a low pass filter can be applied to "
        "node elevation based connected nodes.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2024, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::SwitchArg lowpass_arg(
        "", "lowpass",
        "Applies a lowpass filter to elevation over connected nodes.", false);
    cmd.add(lowpass_arg);
    TCLAP::ValueArg<double> map_static_arg(
        "s", "static",
        "Static elevation to map the input file to. This can be combined with "
        "mapping based on rasters or other meshes to deal with locations where "
        "no corresponding data exists.",
        false, 0, "number");
    cmd.add(map_static_arg);
    TCLAP::ValueArg<double> max_dist_arg(
        "d", "distance",
        "Maximum distance to search for mesh nodes if there is no "
        "corresponding data for input mesh nodes on the mesh it should be "
        "mapped on. (Default value: 1)",
        false, 1, "number");
    cmd.add(max_dist_arg);
    TCLAP::ValueArg<std::string> map_mesh_arg(
        "m", "mesh", "2D mesh file (*.vtu) to map the input file on.", false,
        "", "string");
    cmd.add(map_mesh_arg);
    TCLAP::ValueArg<std::string> map_raster_arg(
        "r", "raster",
        "Raster file (*.asc *.grd *.xyz) to map the input file on.", false, "",
        "string");
    cmd.add(map_raster_arg);
    TCLAP::ValueArg<std::string> output_arg(
        "o", "output", "Output mesh file (*.vtu)", true, "", "string");
    cmd.add(output_arg);
    TCLAP::ValueArg<std::string> input_arg(
        "i", "input", "Input mesh file (*.vtu, *.msh)", true, "", "string");
    cmd.add(input_arg);
    cmd.parse(argc, argv);

    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::IO::readMeshFromFile(input_arg.getValue()));
    if (mesh == nullptr)
    {
        ERR("Error reading mesh file.");
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

    if (!(map_static_arg.isSet() || map_raster_arg.isSet() ||
          map_mesh_arg.isSet()))
    {
        ERR("Nothing to do. Please choose mapping based on a raster or mesh "
            "file, or to a static value.");
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

    if (map_raster_arg.isSet() && map_mesh_arg.isSet())
    {
        ERR("Please select mapping based on *either* a mesh or a raster file.");
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

    // Maps the elevation of mesh nodes to static value
    if (map_static_arg.isSet())
    {
        MeshToolsLib::MeshLayerMapper::mapToStaticValue(
            *mesh, map_static_arg.getValue());
    }

    // Maps the elevation of mesh nodes according to raster
    if (map_raster_arg.isSet())
    {
        std::string const raster_path = map_raster_arg.getValue();
        std::ifstream file_stream(raster_path, std::ifstream::in);
        if (!file_stream.good())
        {
            ERR("Opening raster file {} failed.", raster_path);
#ifdef USE_PETSC
            MPI_Finalize();
#endif
            return EXIT_FAILURE;
        }
        file_stream.close();

        std::unique_ptr<GeoLib::Raster> const raster(
            FileIO::AsciiRasterInterface::readRaster(raster_path));
        MeshToolsLib::MeshLayerMapper::layerMapping(*mesh, *raster, 0, true);
    }

    // Maps the elevation of mesh nodes according to a ground truth mesh
    if (map_mesh_arg.isSet())
    {
        std::unique_ptr<MeshLib::Mesh> ground_truth(
            MeshLib::IO::readMeshFromFile(map_mesh_arg.getValue()));
        if (ground_truth == nullptr)
        {
            ERR("Error reading mesh file.");
#ifdef USE_PETSC
            MPI_Finalize();
#endif
            return EXIT_FAILURE;
        }

        std::vector<MeshLib::Node*> const& nodes = mesh->getNodes();
        MeshLib::MeshElementGrid const grid(*ground_truth);
        auto const edgeLengths = minMaxEdgeLength(mesh->getElements());
        double const max_edge = edgeLengths.second;
        double const max_dist(pow(max_dist_arg.getValue(), 2));

        for (MeshLib::Node* node : nodes)
        {
            MathLib::Point3d min_vol{{(*node)[0] - max_edge,
                                      (*node)[1] - max_edge,
                                      -std::numeric_limits<double>::max()}};
            MathLib::Point3d max_vol{{(*node)[0] + max_edge,
                                      (*node)[1] + max_edge,
                                      std::numeric_limits<double>::max()}};
            std::vector<const MeshLib::Element*> const& elems =
                grid.getElementsInVolume(min_vol, max_vol);
            auto const* element =
                MeshToolsLib::ProjectPointOnMesh::getProjectedElement(elems,
                                                                      *node);
            if (element != nullptr)
            {
                (*node)[2] = MeshToolsLib::ProjectPointOnMesh::getElevation(
                    *element, *node);
            }
            else
            {
                (*node)[2] = getClosestPointElevation(
                    *node, ground_truth->getNodes(), max_dist);
            }
        }
    }

    // a simple lowpass filter for the elevation of mesh nodes using the
    // elevation of each node weighted by 2 and the elevation of each connected
    // node weighted by 1
    if (lowpass_arg.isSet())
    {
        INFO("lowpass");
        const std::size_t nNodes(mesh->getNumberOfNodes());
        std::vector<MeshLib::Node*> nodes(mesh->getNodes());

        std::vector<double> elevation(nNodes);
        for (std::size_t i = 0; i < nNodes; i++)
        {
            elevation[i] = (*nodes[i])[2];
        }

        auto const& connections =
            MeshLib::calculateNodesConnectedByElements(*mesh);
        for (std::size_t i = 0; i < nNodes; i++)
        {
            auto const& conn_nodes(connections[nodes[i]->getID()]);
            const unsigned nConnNodes(conn_nodes.size());
            elevation[i] = (2 * (*nodes[i])[2]);
            for (std::size_t j = 0; j < nConnNodes; ++j)
            {
                elevation[i] += (*conn_nodes[j])[2];
            }
            elevation[i] /= (nConnNodes + 2);
        }

        for (std::size_t i = 0; i < nNodes; i++)
        {
            (*nodes[i])[2] = elevation[i];
        }
    }

    if (MeshLib::IO::writeMeshToFile(*mesh, output_arg.getValue()) != 0)
    {
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

    INFO("Result successfully written.");
#ifdef USE_PETSC
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
