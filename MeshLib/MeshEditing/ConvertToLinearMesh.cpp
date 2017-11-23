/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConvertToLinearMesh.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Hex.h"
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
                  std::vector<MeshLib::Node*> const& vec_new_nodes)
{
    auto const n_base_nodes = T_ELEMENT::n_base_nodes;
    auto** nodes = new MeshLib::Node*[n_base_nodes];
    for (unsigned i=0; i<e->getNumberOfBaseNodes(); i++)
        nodes[i] = const_cast<MeshLib::Node*>(vec_new_nodes[e->getNode(i)->getID()]);
    return new T_ELEMENT(nodes);
}

} // unnamed namespace


std::unique_ptr<MeshLib::Mesh> convertToLinearMesh(MeshLib::Mesh const& org_mesh, std::string const& new_mesh_name)
{
    std::vector<MeshLib::Node*> vec_new_nodes = MeshLib::copyNodeVector(MeshLib::getBaseNodes(org_mesh.getElements()));

    // create new elements with the quadratic nodes
    std::vector<MeshLib::Element*> vec_new_eles;
    for (MeshLib::Element const* e : org_mesh.getElements())
    {
        if (e->getCellType() == MeshLib::CellType::LINE3)
        {
            vec_new_eles.push_back(createLinearElement<MeshLib::Line>(
                e, vec_new_nodes));
        }
        else if (e->getCellType() == MeshLib::CellType::QUAD8)
        {
            vec_new_eles.push_back(createLinearElement<MeshLib::Quad>(
                e, vec_new_nodes));
        }
        else if (e->getCellType() == MeshLib::CellType::HEX20)
        {
            vec_new_eles.push_back(createLinearElement<MeshLib::Hex>(
                e, vec_new_nodes));
        }
        else
        {
            OGS_FATAL("Mesh element type %s is not supported", MeshLib::CellType2String(e->getCellType()).c_str());
        }
    }

    auto new_mesh = std::make_unique<MeshLib::Mesh>(
        new_mesh_name, vec_new_nodes, vec_new_eles,
        org_mesh.getProperties().excludeCopyProperties(
            std::vector<MeshLib::MeshItemType>(1,
                                               MeshLib::MeshItemType::Node)));

    // copy property vectors for nodes
    MeshLib::Properties const& src_properties = org_mesh.getProperties();
    for (auto name : src_properties.getPropertyVectorNames())
    {
        if (!src_properties.existsPropertyVector<double>(name))
            continue;
        auto const* src_prop = src_properties.getPropertyVector<double>(name);
        if (src_prop->getMeshItemType() != MeshLib::MeshItemType::Node)
            continue;

        auto const n_src_comp = src_prop->getNumberOfComponents();
        auto new_prop =
            new_mesh->getProperties().createNewPropertyVector<double>(
                name, MeshLib::MeshItemType::Node, n_src_comp);
        new_prop->resize(new_mesh->getNumberOfNodes() * n_src_comp);

        // copy only base node values
        for (unsigned i=0; i<org_mesh.getNumberOfBaseNodes(); i++)
        {
            for (int j = 0; j < n_src_comp; j++)
                (*new_prop)[i*n_src_comp+j] = (*src_prop)[i*n_src_comp+j];
        }
    }

    return new_mesh;
}

} // end namespace MeshLib

