/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConvertToLinearMesh.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Hex.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Tet.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Elements/Utils.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Properties.h"
#include "MeshLib/PropertyVector.h"
#include "MeshLib/Utils/DuplicateMeshComponents.h"

namespace MeshToolsLib
{
namespace
{
template <typename T_ELEMENT>
T_ELEMENT* createLinearElement(MeshLib::Element const* e,
                               std::vector<MeshLib::Node*> const& vec_new_nodes,
                               std::vector<std::size_t> const& map)
{
    auto const n_base_nodes = T_ELEMENT::n_base_nodes;
    auto** nodes = new MeshLib::Node*[n_base_nodes];
    auto* const* element_nodes = e->getNodes();
    for (unsigned i = 0; i < n_base_nodes; i++)
    {
        auto const n = vec_new_nodes[map[element_nodes[i]->getID()]];
        nodes[i] = const_cast<MeshLib::Node*>(n);
    }
    return new T_ELEMENT(nodes);
}
}  // unnamed namespace

std::unique_ptr<MeshLib::Mesh> convertToLinearMesh(
    MeshLib::Mesh const& mesh, std::string const& new_mesh_name)
{
    auto const& org_elements = mesh.getElements();

    // mark base nodes
    std::vector<bool> marked_base_nodes(mesh.getNodes().size(), false);
    for (auto const org_element : org_elements)
    {
        for (std::size_t k = 0; k < org_element->getNumberOfBaseNodes(); ++k)
        {
            auto const& base_node = *org_element->getNode(k);
            marked_base_nodes[base_node.getID()] = true;
        }
    }

    // construct map and fill new_mesh_nodes
    std::vector<MeshLib::Node*> new_mesh_nodes{static_cast<std::size_t>(
        std::count(begin(marked_base_nodes), end(marked_base_nodes), true))};
    std::size_t base_node_cnt = 0;
    auto const& org_nodes = mesh.getNodes();
    std::vector<std::size_t> base_node_map(org_nodes.size(), -1);
    for (std::size_t k = 0; k < org_nodes.size(); ++k)
    {
        if (marked_base_nodes[k])
        {
            new_mesh_nodes[base_node_cnt] =
                new MeshLib::Node(org_nodes[k]->data(), base_node_cnt);
            base_node_map[k] = base_node_cnt;
            base_node_cnt++;
        }
    }

    // create new elements with the quadratic nodes
    std::vector<MeshLib::Element*> vec_new_eles;
    for (MeshLib::Element const* e : mesh.getElements())
    {
        if (e->getCellType() == MeshLib::CellType::LINE3)
        {
            vec_new_eles.push_back(createLinearElement<MeshLib::Line>(
                e, new_mesh_nodes, base_node_map));
        }
        else if (e->getCellType() == MeshLib::CellType::QUAD8)
        {
            vec_new_eles.push_back(createLinearElement<MeshLib::Quad>(
                e, new_mesh_nodes, base_node_map));
        }
        else if (e->getCellType() == MeshLib::CellType::TRI6)
        {
            vec_new_eles.push_back(createLinearElement<MeshLib::Tri>(
                e, new_mesh_nodes, base_node_map));
        }
        else if (e->getCellType() == MeshLib::CellType::HEX20)
        {
            vec_new_eles.push_back(createLinearElement<MeshLib::Hex>(
                e, new_mesh_nodes, base_node_map));
        }
        else if (e->getCellType() == MeshLib::CellType::TET10)
        {
            vec_new_eles.push_back(createLinearElement<MeshLib::Tet>(
                e, new_mesh_nodes, base_node_map));
        }
        else
        {
            OGS_FATAL("Mesh element type {:s} is not supported",
                      MeshLib::CellType2String(e->getCellType()));
        }
    }

    auto new_mesh = std::make_unique<MeshLib::Mesh>(
        new_mesh_name, new_mesh_nodes, vec_new_eles,
        true /* compute_element_neighbors */,
        mesh.getProperties().excludeCopyProperties(
            std::vector<MeshLib::MeshItemType>(1,
                                               MeshLib::MeshItemType::Node)));

    // copy property vectors for nodes
    for (auto [name, property] : mesh.getProperties())
    {
        if (property->getMeshItemType() != MeshLib::MeshItemType::Node)
        {
            continue;
        }
        auto double_property =
            dynamic_cast<MeshLib::PropertyVector<double>*>(property);
        if (double_property == nullptr)
        {
            continue;
        }
        auto const n_src_comp = double_property->getNumberOfGlobalComponents();
        auto new_prop =
            new_mesh->getProperties().createNewPropertyVector<double>(
                name, MeshLib::MeshItemType::Node, n_src_comp);
        new_prop->resize(new_mesh->getNumberOfNodes() * n_src_comp);

        for (std::size_t k = 0; k < org_nodes.size(); ++k)
        {
            if (!marked_base_nodes[k])
            {
                continue;
            }
            std::copy_n(double_property->begin() + k * n_src_comp,
                        n_src_comp,
                        new_prop->begin() + base_node_map[k] * n_src_comp);
        }
    }
    return new_mesh;
}

}  // namespace MeshToolsLib
