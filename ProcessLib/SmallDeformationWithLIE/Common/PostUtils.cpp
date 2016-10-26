/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "PostUtils.h"

#include <map>
#include <vector>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"

namespace ProcessLib
{
namespace SmallDeformationWithLIE
{

namespace
{

template <typename T>
inline void sort_unique(std::vector<T> &vec)
{
    std::sort(vec.begin(), vec.end());
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}

} // no named namespace

PostProcessTool::PostProcessTool(
    MeshLib::Mesh const& org_mesh,
    std::vector<std::vector<MeshLib::Node*>> const& vec_vec_fracture_nodes,
    std::vector<std::vector<MeshLib::Element*>> const& vec_vec_fracutre_matrix_elements)
    :_org_mesh(org_mesh)
{
    if (!org_mesh.getProperties().hasPropertyVector("displacement")
        || !org_mesh.getProperties().hasPropertyVector("displacement_jump1")
        || !org_mesh.getProperties().hasPropertyVector("levelset1")
        )
    {
        OGS_FATAL("The given mesh does not have relevant properties");
    }

    // clone nodes and elements
    std::vector<MeshLib::Node*> new_nodes(MeshLib::copyNodeVector(org_mesh.getNodes()));
    std::vector<MeshLib::Element*> new_eles(
        MeshLib::copyElementVector(org_mesh.getElements(), new_nodes));

    // duplicate fracture nodes
    for (auto const& vec_fracture_nodes : vec_vec_fracture_nodes)
    {
        for (auto const* org_node : vec_fracture_nodes)
        {
            auto duplicated_node = new MeshLib::Node(org_node->getCoords(), new_nodes.size());
            new_nodes.push_back(duplicated_node);
            assert(_map_dup_newNodeIDs.count(org_node->getID())==0);
            _map_dup_newNodeIDs[org_node->getID()] = duplicated_node->getID();
        }
    }

    // split elements using the new duplicated nodes
    for (unsigned fracture_id=0; fracture_id<=vec_vec_fracutre_matrix_elements.size(); fracture_id++)
    {
        auto const& vec_fracutre_matrix_elements = vec_vec_fracutre_matrix_elements[fracture_id];
        auto const& vec_fracture_nodes = vec_vec_fracture_nodes[fracture_id];
        auto prop_levelset = org_mesh.getProperties().getPropertyVector<double>(
            "levelset" + std::to_string(fracture_id + 1));
        for (auto const* org_e : vec_fracutre_matrix_elements)
        {
            // only matrix elements
            if (org_e->getDimension() != org_mesh.getDimension())
                continue;

            auto const eid = org_e->getID();

            // keep original if the element has levelset=0
            if ((*prop_levelset)[eid] == 0)
                continue;

            // replace fracture nodes with duplicated ones
            MeshLib::Element* e = new_eles[eid];
            for (unsigned i=0; i<e->getNumberOfNodes(); i++)
            {
                // only fracture nodes
                auto itr = _map_dup_newNodeIDs.find(e->getNodeIndex(i));
                if (itr == _map_dup_newNodeIDs.end())
                    continue;

                // check if a node belongs to the particular fracture group
                auto itr2 = std::find_if(
                    vec_fracture_nodes.begin(), vec_fracture_nodes.end(),
                    [&](MeshLib::Node const* node) {
                        return node->getID() == e->getNodeIndex(i);
                    });
                if (itr2 == vec_fracture_nodes.end())
                    continue;

                e->setNode(i, new_nodes[itr->second]);
            }
        }
    }

    // new mesh
    _output_mesh.reset(new MeshLib::Mesh(org_mesh.getName(), new_nodes, new_eles));
    createProperties<int>();
    createProperties<double>();
    copyProperties<int>();
    copyProperties<double>();
    calculateTotalDisplacement(vec_vec_fracture_nodes.size());
}

template <typename T>
void PostProcessTool::createProperties()
{
    MeshLib::Properties const& src_properties = _org_mesh.getProperties();
    for (auto name : src_properties.getPropertyVectorNames())
    {
        auto const* src_prop = src_properties.getPropertyVector<T>(name);
        if (!src_prop)
            continue;
        auto const n_src_comp = src_prop->getNumberOfComponents();
        // convert 2D vector to 3D. Otherwise Paraview Calculator filter does not recognize
        // it as a vector
        auto const n_dest_comp = (n_src_comp==2) ? 3 : n_src_comp;

        if (src_prop->getMeshItemType() == MeshLib::MeshItemType::Node)
        {
            auto new_prop =
                _output_mesh->getProperties().createNewPropertyVector<T>(
                    name, MeshLib::MeshItemType::Node, n_dest_comp);
            new_prop->resize(_output_mesh->getNumberOfNodes() * n_dest_comp);
        }
        else if (src_prop->getMeshItemType() == MeshLib::MeshItemType::Cell)
        {
            auto new_prop =
                _output_mesh->getProperties().createNewPropertyVector<T>(
                    name, MeshLib::MeshItemType::Cell, n_dest_comp);
            new_prop->resize(src_prop->size());
        }
        else
        {
            WARN(
                "Property '%s' cannot be created because its mesh item type is "
                "not supported.",
                name.c_str());
        }
    }
}

template <typename T>
void PostProcessTool::copyProperties()
{
    MeshLib::Properties const& src_properties = _org_mesh.getProperties();
    for (auto name : src_properties.getPropertyVectorNames())
    {
        auto const* src_prop = src_properties.getPropertyVector<T>(name);
        if (!src_prop)
            continue;
        auto* dest_prop = _output_mesh->getProperties().getPropertyVector<T>(name);
        assert(dest_prop);

        if (src_prop->getMeshItemType() == MeshLib::MeshItemType::Node)
        {
            auto const n_src_comp = src_prop->getNumberOfComponents();
            auto const n_dest_comp = dest_prop->getNumberOfComponents();
            // copy existing
            for (unsigned i=0; i<_org_mesh.getNumberOfNodes(); i++)
            {
                for (unsigned j=0; j<n_src_comp; j++)
                    (*dest_prop)[i*n_dest_comp+j] = (*src_prop)[i*n_src_comp+j];
                // set zero for components not existing in the original
                for (unsigned j=n_src_comp; j<n_dest_comp; j++)
                    (*dest_prop)[i*n_dest_comp+j] = 0;
            }
            // copy duplicated
            for (auto itr : _map_dup_newNodeIDs)
            {
                for (unsigned j=0; j<n_dest_comp; j++)
                    (*dest_prop)[itr.second * n_dest_comp + j] =
                        (*dest_prop)[itr.first * n_dest_comp + j];
            }
        }
        else if (src_prop->getMeshItemType() == MeshLib::MeshItemType::Cell)
        {
            std::copy(src_prop->begin(), src_prop->end(), dest_prop->begin());
        }
        else
        {
            WARN(
                "Property '%s' cannot be created because its mesh item type is "
                "not supported.",
                name.c_str());
        }
    }
}

void PostProcessTool::calculateTotalDisplacement(unsigned const n_fractures)
{
    auto const& u =  *_output_mesh->getProperties().getPropertyVector<double>("displacement");
    auto const n_u_comp = u.getNumberOfComponents();
    assert(u.size() == _output_mesh->getNodes().size() * 3);
    auto& total_u =
        *_output_mesh->getProperties().createNewPropertyVector<double>(
            "u", MeshLib::MeshItemType::Node, n_u_comp);
    total_u.resize(u.size());
    for (unsigned i=0; i<_output_mesh->getNodes().size(); i++)
    {
        for (unsigned j=0; j<n_u_comp; j++)
            total_u[i*n_u_comp+j] = u[i*n_u_comp+j];
    }

    for (unsigned fracture_id=0; fracture_id<n_fractures; fracture_id++)
    {
        // nodal value of levelset
        std::vector<double> nodal_levelset(_output_mesh->getNodes().size(), 0.0);
        auto const& ele_levelset =
            *_output_mesh->getProperties().getPropertyVector<double>(
                "levelset" + std::to_string(fracture_id + 1));
        for (MeshLib::Element const* e : _output_mesh->getElements())
        {
            if (e->getDimension() != _output_mesh->getDimension())
                continue;
            const double e_levelset = ele_levelset[e->getID()];
            if (e_levelset == 0)
                continue;

            for (unsigned i=0; i<e->getNumberOfNodes(); i++)
                nodal_levelset[e->getNodeIndex(i)] = e_levelset;
        }

        // update total displacements
        auto const& g =
            *_output_mesh->getProperties().getPropertyVector<double>(
                "displacement_jump" + std::to_string(fracture_id + 1));
        for (unsigned i = 0; i < _output_mesh->getNodes().size(); i++)
        {
            for (unsigned j=0; j<n_u_comp; j++)
                total_u[i*n_u_comp+j] += nodal_levelset[i] * g[i*n_u_comp+j];
        }
    }

}

}  // namespace SmallDeformationWithLIE
}  // namespace ProcessLib
