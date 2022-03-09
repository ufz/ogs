/**
 * \file
 * 2013/13/06 KR Initial implementation
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include <algorithm>
#include <array>
#include <memory>
#include <vector>

#include "BaseLib/Algorithm.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

/**
 * \brief Reverses order of nodes. In particular, this fixes issues between OGS5
 * and OGS6 meshes.
 *
 * \param elements  Mesh elements whose nodes should be reordered
 * \param forced    If true, nodes are reordered for all
 * elements, if false it is first checked if the node order is correct according
 * to OGS6 element definitions.
 */
void reverseNodeOrder(std::vector<MeshLib::Element*>& elements,
                      bool const forced)
{
    std::size_t n_corrected_elements = 0;
    std::size_t const nElements(elements.size());
    for (std::size_t i = 0; i < nElements; ++i)
    {
        if (!forced && elements[i]->testElementNodeOrder())
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

/// Fixes inconsistencies between VTK's and OGS' node order for prism elements.
/// In particular, this fixes issues between OGS6 meshes with and without
/// InSitu-Lib
void fixVtkInconsistencies(std::vector<MeshLib::Element*>& elements)
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

/// Orders the base nodes of each elements before its non-linear nodes.
void reorderNonlinearNodes(MeshLib::Mesh& mesh)
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
                              [](MeshLib::Node* a, MeshLib::Node* b)
                              { return a->getID() < b->getID(); });
    BaseLib::makeVectorUnique(nonlinear_nodes,
                              [](MeshLib::Node* a, MeshLib::Node* b)
                              { return a->getID() < b->getID(); });

    std::vector<MeshLib::Node*>& allnodes =
        const_cast<std::vector<MeshLib::Node*>&>(mesh.getNodes());
    allnodes.clear();

    allnodes.insert(allnodes.end(), base_nodes.begin(), base_nodes.end());
    allnodes.insert(allnodes.end(), nonlinear_nodes.begin(),
                    nonlinear_nodes.end());

    mesh.resetNodeIDs();
}

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Reorders mesh nodes in elements to make old or incorrectly ordered "
        "meshes compatible with OGS6.\n"
        "Three options are available:\n"
        "Method 0: Reversing order of nodes for all elements.\n"
        "Method 1: Reversing order of nodes unless it's perceived correct by "
        "OGS6 standards. This is the default selection.\n"
        "Method 2: Fixing node ordering issues between VTK and OGS6 (only "
        "applies to prism-elements)\n"
        "Method 3: Re-ordering of mesh node vector such that all base nodes "
        "are sorted before all nonlinear nodes.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2022, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    std::vector<int> method_ids{0, 1, 2, 3};
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
    if (!method_arg.isSet() || method_arg.getValue() < 2)
    {
        bool const forced = (method_arg.getValue() == 0);
        reverseNodeOrder(
            const_cast<std::vector<MeshLib::Element*>&>(mesh->getElements()),
            forced);
    }
    else if (method_arg.getValue() == 2)
    {
        fixVtkInconsistencies(
            const_cast<std::vector<MeshLib::Element*>&>(mesh->getElements()));
    }
    else if (method_arg.getValue() == 3)
    {
        reorderNonlinearNodes(*mesh);
    }

    MeshLib::IO::writeMeshToFile(*mesh, output_mesh_arg.getValue());

    INFO("VTU file written.");

    return EXIT_SUCCESS;
}
