/**
 * \file
 * 2012/03/07 KR Initial implementation
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <memory>
#include <string>

#ifdef USE_PETSC
#include <mpi.h>
#endif

#include "BaseLib/FileTools.h"
#include "GeoLib/AABB.h"
#include "InfoLib/GitInfo.h"
#include "MathLib/MathTools.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/ProjectPointOnMesh.h"
#include "MeshLib/MeshSearch/MeshElementGrid.h"
#include "MeshLib/Node.h"

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
    std::vector<std::string> keywords;
    keywords.emplace_back("-ALL");
    keywords.emplace_back("-MESH");
    keywords.emplace_back("-LOWPASS");

    if (argc < 3)
    {
        INFO(
            "Moves mesh nodes and connected elements either by a given value "
            "or based on a list.\n");
        INFO("Usage: {:s} <msh-file.msh> <keyword> [<value1>] [<value2>]",
             argv[0]);
        INFO("Available keywords:");
        INFO(
            "\t-MESH <value1> <value2> : changes the elevation of mesh nodes "
            "based on a second mesh <value1> with a search range of <value2>.");
        INFO(
            "\t-LOWPASS : applies a lowpass filter over node elevation using "
            "directly connected nodes.");
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

    const std::string msh_name(argv[1]);
    const std::string current_key(argv[2]);
    std::string const ext(BaseLib::getFileExtension(msh_name));
    if (!(ext == ".msh" || ext == ".vtu"))
    {
        ERR("Error: Parameter 1 must be a mesh-file (*.msh / *.vtu).");
        INFO("Usage: {:s} <msh-file.gml> <keyword> <value>", argv[0]);
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

    bool const is_keyword = std::any_of(keywords.begin(), keywords.end(),
                                        [current_key](auto const& keyword)
                                        { return current_key == keyword; });

    if (!is_keyword)
    {
        ERR("Keyword not recognised. Available keywords:");
        for (auto const& keyword : keywords)
            INFO("\t{:s}", keyword);
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::IO::readMeshFromFile(msh_name));
    if (mesh == nullptr)
    {
        ERR("Error reading mesh file.");
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return 1;
    }

    // Start keyword-specific selection of nodes
    // maps the elevation of mesh nodes according to a ground truth mesh
    // whenever nodes exist within max_dist
    if (current_key == "-MESH")
    {
        if (argc < 4)
        {
            ERR("Missing parameter...");
#ifdef USE_PETSC
            MPI_Finalize();
#endif
            return EXIT_FAILURE;
        }
        const std::string value(argv[3]);
        double max_dist(pow(strtod(argv[4], 0), 2));
        std::unique_ptr<MeshLib::Mesh> ground_truth(
            MeshLib::IO::readMeshFromFile(value));
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
        double const max_edge(mesh->getMaxEdgeLength());

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
                MeshLib::ProjectPointOnMesh::getProjectedElement(elems, *node);
            (*node)[2] =
                (element != nullptr)
                    ? MeshLib::ProjectPointOnMesh::getElevation(*element, *node)
                    : getClosestPointElevation(*node, ground_truth->getNodes(),
                                               max_dist);
        }
    }

    // a simple lowpass filter for the elevation of mesh nodes using the
    // elevation of each node weighted by 2 and the elevation of each connected
    // node weighted by 1
    if (current_key == "-LOWPASS")
    {
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
    /**** add other keywords here ****/

    std::string const new_mesh_name(msh_name.substr(0, msh_name.length() - 4) +
                                    "_new.vtu");
    if (MeshLib::IO::writeMeshToFile(*mesh, new_mesh_name) != 0)
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
