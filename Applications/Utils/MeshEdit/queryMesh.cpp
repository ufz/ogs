/**
 * \file
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

#include <array>
#include <memory>
#include <sstream>
#include <string>

#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Query mesh information.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2024, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::UnlabeledValueArg<std::string> mesh_arg(
        "mesh-file", "input mesh file", true, "", "string");
    cmd.add(mesh_arg);
    TCLAP::MultiArg<std::size_t> eleId_arg("e", "element-id", "element ID",
                                           false, "number");
    cmd.add(eleId_arg);
    TCLAP::MultiArg<std::size_t> nodeId_arg("n", "node-id", "node ID", false,
                                            "number");
    cmd.add(nodeId_arg);
    TCLAP::SwitchArg showNodeWithMaxEle_arg(
        "", "show-node-with-max-elements",
        "show a node having the max number of connected elements", false);
    cmd.add(showNodeWithMaxEle_arg);

    cmd.parse(argc, argv);

#ifdef USE_PETSC
    MPI_Init(&argc, &argv);
#endif

    const std::string filename(mesh_arg.getValue());

    // read the mesh file
    auto const mesh =
        std::unique_ptr<MeshLib::Mesh>(MeshLib::IO::readMeshFromFile(
            filename, true /* compute_element_neighbors */));
    if (!mesh)
    {
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

    std::vector<std::size_t> selected_node_ids;
    if (showNodeWithMaxEle_arg.getValue())
    {
        auto itr = std::max_element(
            mesh->getNodes().begin(), mesh->getNodes().end(),
            [&mesh](MeshLib::Node* i, MeshLib::Node* j)
            {
                return mesh->getElementsConnectedToNode(*i).size() <
                       mesh->getElementsConnectedToNode(*j).size();
            });
        if (itr != mesh->getNodes().end())
        {
            MeshLib::Node* node = *itr;
            selected_node_ids.push_back(node->getID());
        }
    }
    selected_node_ids.insert(selected_node_ids.end(),
                             nodeId_arg.getValue().begin(),
                             nodeId_arg.getValue().end());

    auto const materialIds = materialIDs(*mesh);
    for (auto ele_id : eleId_arg.getValue())
    {
        std::stringstream out;
        out << std::scientific
            << std::setprecision(std::numeric_limits<double>::digits10);
        out << "--------------------------------------------------------"
            << std::endl;
        auto* ele = mesh->getElement(ele_id);
        out << "# Element " << ele->getID() << std::endl;
        out << "Type : " << CellType2String(ele->getCellType()) << std::endl;
        if (materialIds)
        {
            out << "Mat ID : " << (*materialIds)[ele_id] << std::endl;
        }
        out << "Nodes: " << std::endl;
        for (unsigned i = 0; i < ele->getNumberOfNodes(); i++)
        {
            out << ele->getNode(i)->getID() << " " << *ele->getNode(i)
                << std::endl;
        }
        out << "Content: " << ele->getContent() << std::endl;
        out << "Neighbors: ";
        for (unsigned i = 0; i < ele->getNumberOfNeighbors(); i++)
        {
            if (ele->getNeighbor(i))
            {
                out << ele->getNeighbor(i)->getID() << " ";
            }
            else
            {
                out << "none ";
            }
        }
        out << std::endl;
        INFO("{:s}", out.str());
    }

    auto const& connections = MeshLib::calculateNodesConnectedByElements(*mesh);
    for (auto node_id : selected_node_ids)
    {
        std::stringstream out;
        out << std::scientific
            << std::setprecision(std::numeric_limits<double>::digits10);
        out << "--------------------------------------------------------"
            << std::endl;
        MeshLib::Node const* node = mesh->getNode(node_id);
        out << "# Node " << node->getID() << std::endl;
        out << "Coordinates: " << *node << std::endl;
        out << "Connected elements ("
            << mesh->getElementsConnectedToNode(*node).size() << "): ";
        for (auto ele : mesh->getElementsConnectedToNode(*node))
        {
            out << ele->getID() << " ";
        }
        out << std::endl;
        out << "Connected nodes (" << connections[node->getID()].size()
            << "): ";
        for (auto nd : connections[node->getID()])
        {
            out << nd->getID() << " ";
        }
        out << std::endl;
        INFO("{:s}", out.str());
    }
#ifdef USE_PETSC
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
