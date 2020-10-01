/**
 * \file
 * 2013/13/06 KR Initial implementation
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <array>
#include <algorithm>
#include <memory>
#include <vector>

#include <tclap/CmdLine.h>

#include "BaseLib/Algorithm.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

/// Re-ordering mesh elements to correct Data Explorer 5 meshes to work with Data Explorer 6.
void reorderNodes(std::vector<MeshLib::Element*>& elements)
{
    std::size_t n_corrected_elements = 0;
    std::size_t const nElements(elements.size());
    for (std::size_t i = 0; i < nElements; ++i)
    {
        if (elements[i]->testElementNodeOrder())
        {
            continue;
        }
        n_corrected_elements++;

        const unsigned nElemNodes(elements[i]->getNumberOfBaseNodes());
        std::vector<MeshLib::Node*> nodes(elements[i]->getNodes(),
                                          elements[i]->getNodes() + nElemNodes);

        switch (elements[i]->getGeomType())
        {
            case MeshLib::MeshElemType::TETRAHEDRON:
                for (std::size_t j = 0; j < 4; ++j)
                {
                    elements[i]->setNode(j, nodes[(j + 1) % 4]);
                }
                break;
            case MeshLib::MeshElemType::PYRAMID:
                elements[i]->setNode(0, nodes[1]);
                elements[i]->setNode(1, nodes[0]);
                elements[i]->setNode(2, nodes[3]);
                elements[i]->setNode(3, nodes[2]);
                break;
            case MeshLib::MeshElemType::PRISM:
                for (std::size_t j = 0; j < 3; ++j)
                {
                    elements[i]->setNode(j, nodes[j + 3]);
                    elements[i]->setNode(j + 3, nodes[j]);
                }
                break;
            case MeshLib::MeshElemType::HEXAHEDRON:
                for (std::size_t j = 0; j < 4; ++j)
                {
                    elements[i]->setNode(j, nodes[j + 4]);
                    elements[i]->setNode(j + 4, nodes[j]);
                }
                break;
            default:
                for (std::size_t j = 0; j < nElemNodes; ++j)
                {
                    elements[i]->setNode(j, nodes[nElemNodes - j - 1]);
                }
        }
    }

    INFO("Corrected {:d} elements.", n_corrected_elements);
}

/// Re-ordering prism elements to correct OGS6 meshes with and without InSitu-Lib
void reorderNodes2(std::vector<MeshLib::Element*>& elements)
{
    std::size_t const nElements(elements.size());
    for (std::size_t i = 0; i < nElements; ++i)
    {
        const unsigned nElemNodes(elements[i]->getNumberOfBaseNodes());
        std::vector<MeshLib::Node*> nodes(elements[i]->getNodes(),
                                          elements[i]->getNodes() + nElemNodes);

        for (std::size_t j = 0; j < nElemNodes; ++j)
        {
            if (elements[i]->getGeomType() == MeshLib::MeshElemType::PRISM)
            {
                for (std::size_t k = 0; k < 3; ++k)
                {
                    elements[i]->setNode(k, nodes[k + 3]);
                    elements[i]->setNode(k + 3, nodes[k]);
                }
                break;
            }
        }
    }
}

void reorderNonlinearNodes(MeshLib::Mesh &mesh)
{
    std::vector<MeshLib::Node*> base_nodes;
    std::vector<MeshLib::Node*> nonlinear_nodes;
    for (MeshLib::Element const* e : mesh.getElements())
    {
        for (unsigned i = 0; i < e->getNumberOfBaseNodes(); i++)
        {
            base_nodes.push_back(const_cast<MeshLib::Node*>(e->getNode(i)));
        }
        for (unsigned i = e->getNumberOfBaseNodes(); i < e->getNumberOfNodes();
             i++)
        {
            nonlinear_nodes.push_back(
                const_cast<MeshLib::Node*>(e->getNode(i)));
        }
    }

    BaseLib::makeVectorUnique(base_nodes,
                              [](MeshLib::Node* a, MeshLib::Node* b) {
                                  return a->getID() < b->getID();
                              });
    BaseLib::makeVectorUnique(nonlinear_nodes,
                              [](MeshLib::Node* a, MeshLib::Node* b) {
                                  return a->getID() < b->getID();
                              });

    std::vector<MeshLib::Node*>& allnodes =
        const_cast<std::vector<MeshLib::Node*>&>(mesh.getNodes());
    allnodes.clear();

    allnodes.insert(allnodes.end(), base_nodes.begin(), base_nodes.end());
    allnodes.insert(allnodes.end(), nonlinear_nodes.begin(), nonlinear_nodes.end());

    mesh.resetNodeIDs();
}

int main (int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Reorders mesh nodes in elements to make old or incorrectly ordered "
        "meshes compatible with OGS6.\n"
        "Three options are available:\n"
        "Method 1: Re-ordering between DataExplorer 5 and DataExplorer 6 "
        "(mostly reverses order of the nodes for all element types\n"
        "Method 2: Re-ordering after introducing InSitu-Lib to OGS6 (only "
        "adjusts to top and bottom surfaces of prism elements\n"
        "Method 3: Re-ordering of mesh node vector such that all base nodes "
        "are sorted before all nonlinear nodes.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2020, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    std::vector<int> method_ids{1, 2, 3};
    TCLAP::ValuesConstraint<int> allowed_values(method_ids);
    TCLAP::ValueArg<int> method_arg("m", "method",
                                    "reordering method selection", false, 1,
                                    &allowed_values);
    cmd.add(method_arg);
    TCLAP::ValueArg<std::string> output_mesh_arg(
        "o", "output_mesh", "the name of the output mesh file", true, "",
        "filename");
    cmd.add(output_mesh_arg);
    TCLAP::ValueArg<std::string> input_mesh_arg(
        "i", "input_mesh", "the name of the input mesh file", true, "",
        "filename");
    cmd.add(input_mesh_arg);
    cmd.parse(argc, argv);

    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::IO::readMeshFromFile(input_mesh_arg.getValue()));

    if (!mesh)
    {
        return EXIT_FAILURE;
    }

    INFO("Reordering nodes... ");
    if (!method_arg.isSet() || method_arg.getValue() == 1)
    {
        reorderNodes(const_cast<std::vector<MeshLib::Element*>&>(mesh->getElements()));
    }
    else if (method_arg.getValue() == 2)
    {
        reorderNodes2(const_cast<std::vector<MeshLib::Element*>&>(mesh->getElements()));
    }
    else if (method_arg.getValue() == 3)
    {
        reorderNonlinearNodes(*mesh);
    }

    MeshLib::IO::writeMeshToFile(*mesh, output_mesh_arg.getValue());

    INFO("VTU file written.");

    return EXIT_SUCCESS;
}



