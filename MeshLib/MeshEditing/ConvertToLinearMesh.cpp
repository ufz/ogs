/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
#include "MeshLib/Node.h"
#include "MeshLib/Properties.h"
#include "MeshLib/PropertyVector.h"

namespace MeshLib
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
    MeshLib::Mesh const& org_mesh, std::string const& new_mesh_name)
{
    auto const& org_elements = org_mesh.getElements();
    std::vector<MeshLib::Node*> vec_new_nodes =
        MeshLib::copyNodeVector(MeshLib::getBaseNodes(org_elements));

    // map from old node ids (e->getNode(i)) to new node ids (vec_new_nodes).
    std::vector<std::size_t> map(org_mesh.getNumberOfNodes(), -1);
    for (std::size_t i = 0; i < vec_new_nodes.size(); ++i)
    {
        auto const it = find_if(
            begin(org_mesh.getNodes()), end(org_mesh.getNodes()),
            [node_i = vec_new_nodes[i]](Node* const org_node)
            {
                return *node_i ==
                       *org_node;  // coordinate comparison up to epsilon
            });
        if (it == end(org_mesh.getNodes()))
        {
            OGS_FATAL("A base node");
        }
        map[(*it)->getID()] = i;
    }



    // create new elements with the quadratic nodes
    std::vector<MeshLib::Element*> vec_new_eles;
    for (MeshLib::Element const* e : org_mesh.getElements())
    {
        if (e->getCellType() == MeshLib::CellType::LINE3)
        {
            vec_new_eles.push_back(
                createLinearElement<MeshLib::Line>(e, vec_new_nodes, map));
        }
        else if (e->getCellType() == MeshLib::CellType::QUAD8)
        {
            vec_new_eles.push_back(
                createLinearElement<MeshLib::Quad>(e, vec_new_nodes, map));
        }
        else if (e->getCellType() == MeshLib::CellType::TRI6)
        {
            vec_new_eles.push_back(
                createLinearElement<MeshLib::Tri>(e, vec_new_nodes, map));
        }
        else if (e->getCellType() == MeshLib::CellType::HEX20)
        {
            vec_new_eles.push_back(
                createLinearElement<MeshLib::Hex>(e, vec_new_nodes, map));
        }
        else if (e->getCellType() == MeshLib::CellType::TET10)
        {
            vec_new_eles.push_back(
                createLinearElement<MeshLib::Tet>(e, vec_new_nodes, map));
        }
        else
        {
            OGS_FATAL("Mesh element type {:s} is not supported",
                      MeshLib::CellType2String(e->getCellType()));
        }
    }

    auto new_mesh = std::make_unique<MeshLib::Mesh>(
        new_mesh_name, vec_new_nodes, vec_new_eles,
        org_mesh.getProperties().excludeCopyProperties(
            std::vector<MeshLib::MeshItemType>(1,
                                               MeshLib::MeshItemType::Node)));

    // copy property vectors for nodes
    for (auto [name, property] : org_mesh.getProperties())
    {
        if (property->getMeshItemType() != MeshLib::MeshItemType::Node)
        {
            continue;
        }
        auto double_property = dynamic_cast<PropertyVector<double>*>(property);
        if (double_property == nullptr)
        {
            continue;
        }
        auto const n_src_comp = double_property->getNumberOfGlobalComponents();
        auto new_prop =
            new_mesh->getProperties().createNewPropertyVector<double>(
                name, MeshLib::MeshItemType::Node, n_src_comp);
        new_prop->resize(new_mesh->getNumberOfNodes() * n_src_comp);

        // copy only base node values
        for (unsigned i = 0; i < org_mesh.getNumberOfBaseNodes(); i++)
        {
            for (int j = 0; j < n_src_comp; j++)
            {
                (*new_prop)[i * n_src_comp + j] =
                    (*double_property)[i * n_src_comp + j];
            }
        }
    }

    return new_mesh;
}

}  // end namespace MeshLib
