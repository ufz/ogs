/**
 * \file moveMeshNodes.cpp
 * 2012/03/07 KR Initial implementation
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <memory>
#include <string>

#include "Applications/ApplicationsLib/LogogSetup.h"
#include "BaseLib/BuildInfo.h"
#include "BaseLib/FileTools.h"
#include "GeoLib/AABB.h"
#include "MathLib/MathTools.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/MeshEditing/ProjectPointOnMesh.h"

bool setClosestPointElevation(MeshLib::Node& p,
                              std::vector<MeshLib::Node*> const& nodes,
                              double const& max_dist)
{
    double sqr_shortest_dist (max_dist);
    bool pnt_found (false);
    for (MeshLib::Node* node : nodes)
    {
        double sqr_dist = (p[0] - (*node)[0]) * (p[0] - (*node)[0]) +
                          (p[1] - (*node)[1]) * (p[1] - (*node)[1]);
        if (sqr_dist < sqr_shortest_dist)
        {
            sqr_shortest_dist = sqr_dist;
            p[2] = (*node)[2];
            pnt_found = true;
        }
    }
    return pnt_found;
}

int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    std::vector<std::string> keywords;
    keywords.push_back("-ALL");
    keywords.push_back("-MESH");
    keywords.push_back("-LOWPASS");

    if (argc < 3)
    {
        INFO(
            "Moves mesh nodes and connected elements either by a given value "
            "or based on a list.\n");
        INFO("Usage: %s <msh-file.msh> <keyword> [<value1>] [<value2>]",
             argv[0]);
        INFO("Available keywords:");
        INFO(
            "\t-ALL <value1> <value2> : changes the elevation of all mesh "
            "nodes by <value2> in direction <value1> [x,y,z].");
        INFO(
            "\t-MESH <value1> <value2> : changes the elevation of mesh nodes "
            "based on a second mesh <value1> with a search range of <value2>.");
        INFO(
            "\t-LOWPASS : applies a lowpass filter over node elevation using "
            "directly connected nodes.");
        return EXIT_FAILURE;
    }

    const std::string msh_name(argv[1]);
    const std::string current_key(argv[2]);
    std::string const ext (BaseLib::getFileExtension(msh_name));
    if (!(ext == "msh" || ext == "vtu"))
    {
        ERR("Error: Parameter 1 must be a mesh-file (*.msh / *.vtu).");
        INFO("Usage: %s <msh-file.gml> <keyword> <value>", argv[0]);
        return EXIT_FAILURE;
    }

    bool is_keyword(false);
    for (auto& keyword : keywords)
    {
        if (current_key == keyword)
        {
            is_keyword = true;
            break;
        }
    }

    if (!is_keyword)
    {
        ERR("Keyword not recognised. Available keywords:");
        for (auto const& keyword : keywords)
            INFO("\t%s", keyword.c_str());
        return EXIT_FAILURE;
    }

    std::unique_ptr<MeshLib::Mesh> mesh (MeshLib::IO::readMeshFromFile(msh_name));
    if (mesh == nullptr)
    {
        ERR ("Error reading mesh file.");
        return 1;
    }

    // Start keyword-specific selection of nodes

    // moves the elevation of all nodes by value
    if (current_key == "-ALL")
    {
        if (argc < 5)
        {
            ERR("Missing parameter...");
            return EXIT_FAILURE;
        }
        const std::string dir(argv[3]);
        unsigned idx = (dir == "x") ? 0 : (dir == "y") ? 1 : 2;
        const double value(strtod(argv[4],0));
        INFO("Moving all mesh nodes by %g in direction %d (%s)...", value, idx,
             dir.c_str());
        //double value(-10);
        const std::size_t nNodes(mesh->getNumberOfNodes());
        std::vector<MeshLib::Node*> nodes (mesh->getNodes());
        for (std::size_t i=0; i<nNodes; i++)
        {
            (*nodes[i])[idx] += value;
        }
    }

    // maps the elevation of mesh nodes according to a ground truth mesh whenever nodes exist within max_dist
    if (current_key == "-MESH")
    {
        if (argc < 4)
        {
            ERR("Missing parameter...");
            return EXIT_FAILURE;
        }
        const std::string value (argv[3]);
        double max_dist(pow(strtod(argv[4],0), 2));
        double offset (0.0); // additional offset for elevation (should be 0)
        std::unique_ptr<MeshLib::Mesh> ground_truth (MeshLib::IO::readMeshFromFile(value));
        if (ground_truth == nullptr)
        {
            ERR ("Error reading mesh file.");
            return EXIT_FAILURE;
        }

        std::vector<MeshLib::Node*> const& nodes = mesh->getNodes();
        MeshLib::ProjectPointOnMesh::project(
            *ground_truth, nodes, std::numeric_limits<double>::max());
        std::size_t pnts_not_found (0);
        for (MeshLib::Node* node : nodes)
        {
            if ((*node)[2] == std::numeric_limits<double>::max())
            {
                if (!setClosestPointElevation(*node, ground_truth->getNodes(), max_dist))
                    pnts_not_found++;
            }
        }
        if (pnts_not_found > 0)
            WARN ("For %d points no corresponding elevation found.", pnts_not_found);
    }

    // a simple lowpass filter for the elevation of mesh nodes using the elevation of each node
    // weighted by 2 and the elevation of each connected node weighted by 1
    if (current_key == "-LOWPASS")
    {
        const std::size_t nNodes(mesh->getNumberOfNodes());
        std::vector<MeshLib::Node*> nodes (mesh->getNodes());

        std::vector<double> elevation(nNodes);
        for (std::size_t i = 0; i < nNodes; i++)
        {
            elevation[i] = (*nodes[i])[2];
        }

        for (std::size_t i=0; i<nNodes; i++)
        {
            const std::vector<MeshLib::Node*> conn_nodes (nodes[i]->getConnectedNodes());
            const unsigned nConnNodes (conn_nodes.size());
            elevation[i] = (2*(*nodes[i])[2]);
            for (std::size_t j = 0; j < nConnNodes; ++j)
            {
                elevation[i] += (*conn_nodes[j])[2];
            }
            elevation[i] /= (nConnNodes+2);
        }

        for (std::size_t i = 0; i < nNodes; i++)
        {
            (*nodes[i])[2] = elevation[i];
        }
    }
    /**** add other keywords here ****/

    std::string const new_mesh_name (msh_name.substr(0, msh_name.length() - 4) + "_new.vtu");
    if (MeshLib::IO::writeMeshToFile(*mesh, new_mesh_name) != 0)
    {
        return EXIT_FAILURE;
    }

    INFO("Result successfully written.");
    return EXIT_SUCCESS;
}
