/**
 * \file moveMeshNodes.cpp
 * 2012/03/07 KR Initial implementation
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <memory>
#include <string>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/FileTools.h"

#include "GeoLib/AABB.h"

#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

#include "MathLib/MathTools.h"

int find_closest_point(MeshLib::Node const*const point, std::vector<MeshLib::Node*> const& nodes, double const& max_dist)
{
    const std::size_t nNodes (nodes.size());
    double sqr_shortest_dist (max_dist*2);
    int idx = (sqr_shortest_dist<max_dist) ? 0 : -1;
    const MeshLib::Node p (*point);

    for (unsigned i=0; i<nNodes; i++)
    {
        double sqr_dist ((p[0]-(*nodes[i])[0])*(p[0]-(*nodes[i])[0]));
        if (sqr_dist < max_dist)
        {
            sqr_dist += ((p[1]-(*nodes[i])[1])*(p[1]-(*nodes[i])[1]));
            if (sqr_dist < max_dist && sqr_dist < sqr_shortest_dist)
            {
                sqr_shortest_dist = sqr_dist;
                idx = i;
            }
        }
    }

    return idx;
}

bool containsPoint(MeshLib::Node const& pnt, MathLib::Point3d const& min,
    MathLib::Point3d const& max)
{
    if (pnt[0] < min[0] || max[0] < pnt[0]) return false;
    if (pnt[1] < min[1] || max[1] < pnt[1]) return false;
    return true;
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
        //for (std::size_t i=0; i<keywords.size(); i++)
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
    for (auto & keyword : keywords)
        if (current_key.compare(keyword)==0)
        {
            is_keyword = true;
            break;
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
    if (current_key.compare("-ALL")==0)
    {
        if (argc < 5)
        {
            ERR("Missing parameter...");
            return EXIT_FAILURE;
        }
        const std::string dir(argv[3]);
        unsigned idx = (dir.compare("x") == 0) ? 0 : (dir.compare("y") == 0) ? 1 : 2;
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
    if (current_key.compare("-MESH")==0)
    {
        if (argc < 5)
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

        const std::vector<MeshLib::Node*>& ground_truth_nodes (ground_truth->getNodes());
        GeoLib::AABB bounding_box(ground_truth_nodes.begin(), ground_truth_nodes.end());
        MathLib::Point3d const& min(bounding_box.getMinPoint());
        MathLib::Point3d const& max(bounding_box.getMaxPoint());

        const std::size_t nNodes(mesh->getNumberOfNodes());
        std::vector<MeshLib::Node*> nodes (mesh->getNodes());

        for (std::size_t i=0; i<nNodes; i++)
        {
            bool is_inside (containsPoint(*nodes[i], min, max));
            if (is_inside)
            {
                int idx = find_closest_point(nodes[i], ground_truth_nodes, max_dist);
                if (idx>=0)
                    (*nodes[i])[2] = (*(ground_truth_nodes[idx]))[2]-offset;
            }
        }
    }

    // a simple lowpass filter for the elevation of mesh nodes using the elevation of each node
    // weighted by 2 and the elevation of each connected node weighted by 1
    if (current_key.compare("-LOWPASS")==0)
    {
        const std::size_t nNodes(mesh->getNumberOfNodes());
        std::vector<MeshLib::Node*> nodes (mesh->getNodes());

        std::vector<double> elevation(nNodes);
        for (std::size_t i=0; i<nNodes; i++)
            elevation[i] = (*nodes[i])[2];

        for (std::size_t i=0; i<nNodes; i++)
        {
            const std::vector<MeshLib::Node*> conn_nodes (nodes[i]->getConnectedNodes());
            const unsigned nConnNodes (conn_nodes.size());
            elevation[i] = (2*(*nodes[i])[2]);
            for (std::size_t j=0; j<nConnNodes; ++j)
                elevation[i] += (*conn_nodes[j])[2];
            elevation[i] /= (nConnNodes+2);
        }

        for (std::size_t i=0; i<nNodes; i++)
            (*nodes[i])[2] = elevation[i];
    }
    /**** add other keywords here ****/

    MeshLib::IO::VtuInterface vtu (mesh.get(), 0, false);
    vtu.writeToFile(msh_name.substr(0, msh_name.length() - 4) + "_new.vtu");
    return EXIT_SUCCESS;
}
