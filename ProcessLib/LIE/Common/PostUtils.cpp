/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "PostUtils.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
#include "MeshLib/Node.h"

namespace ProcessLib
{
namespace LIE
{
namespace
{
bool includesNodeID(std::vector<MeshLib::Node*> const& vec_nodes,
                    std::size_t node_id)
{
    return std::any_of(vec_nodes.begin(), vec_nodes.end(),
                       [&](MeshLib::Node const* node)
                       { return node->getID() == node_id; });
}

std::vector<int> const& getMaterialIdsForNode(
    std::vector<std::pair<std::size_t, std::vector<int>>> const&
        vec_nodeID_matIDs,
    std::size_t nodeID)
{
    auto itr =
        std::find_if(vec_nodeID_matIDs.begin(), vec_nodeID_matIDs.end(),
                     [&](std::pair<std::size_t, std::vector<int>> const& entry)
                     { return entry.first == nodeID; });
    assert(itr != vec_nodeID_matIDs.end());
    return itr->second;
};
}  // namespace

PostProcessTool::PostProcessTool(
    MeshLib::Mesh const& org_mesh,
    std::vector<int> const& vec_fracture_mat_IDs,
    std::vector<std::vector<MeshLib::Node*>> const& vec_vec_fracture_nodes,
    std::vector<std::vector<MeshLib::Element*>> const&
        vec_vec_fracture_matrix_elements,
    std::vector<std::pair<std::size_t, std::vector<int>>> const&
        vec_branch_nodeID_matIDs,
    std::vector<std::pair<std::size_t, std::vector<int>>> const&
        vec_junction_nodeID_matIDs)
    : _org_mesh(org_mesh)
{
    if (!org_mesh.getProperties().hasPropertyVector("displacement") ||
        !org_mesh.getProperties().hasPropertyVector("displacement_jump1") ||
        !org_mesh.getProperties().hasPropertyVector("levelset1"))
    {
        OGS_FATAL("The given mesh does not have relevant properties");
    }

    // clone nodes and elements
    std::vector<MeshLib::Node*> new_nodes(
        MeshLib::copyNodeVector(org_mesh.getNodes()));
    std::vector<MeshLib::Element*> new_eles(
        MeshLib::copyElementVector(org_mesh.getElements(), new_nodes));

    // duplicate fracture nodes (two dup. nodes created at a branch)
    for (auto const& vec_fracture_nodes : vec_vec_fracture_nodes)
    {
        for (auto const* org_node : vec_fracture_nodes)
        {
            auto duplicated_node =
                new MeshLib::Node(org_node->data(), new_nodes.size());
            new_nodes.push_back(duplicated_node);
            _map_dup_newNodeIDs[org_node->getID()].push_back(
                duplicated_node->getID());
        }
    }
    // at a junction, generate one more duplicated node (total 4 nodes)
    for (auto& entry : vec_junction_nodeID_matIDs)
    {
        auto* org_node = org_mesh.getNode(entry.first);
        auto duplicated_node =
            new MeshLib::Node(org_node->data(), new_nodes.size());
        new_nodes.push_back(duplicated_node);
        _map_dup_newNodeIDs[org_node->getID()].push_back(
            duplicated_node->getID());
    }

    // split elements using the new duplicated nodes
    for (unsigned fracture_id = 0;
         fracture_id < vec_vec_fracture_matrix_elements.size();
         fracture_id++)
    {
        auto const frac_matid = vec_fracture_mat_IDs[fracture_id];
        auto const& vec_fracture_matrix_elements =
            vec_vec_fracture_matrix_elements[fracture_id];
        auto const& vec_fracture_nodes = vec_vec_fracture_nodes[fracture_id];
        auto prop_levelset = org_mesh.getProperties().getPropertyVector<double>(
            "levelset" + std::to_string(fracture_id + 1));
        for (auto const* org_e : vec_fracture_matrix_elements)
        {
            // only matrix elements
            if (org_e->getDimension() != org_mesh.getDimension())
            {
                continue;
            }

            auto const eid = org_e->getID();
            // keep original if the element has levelset<=0
            if ((*prop_levelset)[eid] <= 0)
            {
                continue;
            }

            // replace fracture nodes with duplicated ones
            MeshLib::Element* e = new_eles[eid];
            for (unsigned i = 0; i < e->getNumberOfNodes(); i++)
            {
                const auto node_id = getNodeIndex(*e, i);
                if (!includesNodeID(vec_fracture_nodes, node_id))
                {
                    continue;
                }

                // list of duplicated node IDs
                auto itr = _map_dup_newNodeIDs.find(node_id);
                if (itr == _map_dup_newNodeIDs.end())
                {
                    continue;
                }
                const auto& dup_newNodeIDs = itr->second;

                // choose new node id
                unsigned new_node_id = 0;
                if (dup_newNodeIDs.size() == 1)
                {
                    // non-intersected nodes
                    new_node_id = dup_newNodeIDs[0];
                }
                else if (dup_newNodeIDs.size() == 2)
                {
                    // branch nodes
                    const auto& br_matids = getMaterialIdsForNode(
                        vec_branch_nodeID_matIDs, node_id);
                    auto const frac2_matid = BaseLib::findFirstNotEqualElement(
                        br_matids, frac_matid);
                    assert(frac2_matid);
                    auto const frac2_id =
                        BaseLib::findIndex(vec_fracture_mat_IDs, *frac2_matid);
                    assert(frac2_id != std::numeric_limits<std::size_t>::max());
                    auto prop_levelset2 =
                        org_mesh.getProperties().getPropertyVector<double>(
                            "levelset" + std::to_string(frac2_id + 1));
                    std::size_t pos = 0;
                    if ((*prop_levelset2)[eid] <= 0)
                    {
                        // index of this fracture
                        pos = BaseLib::findIndex(br_matids, frac_matid);
                    }
                    else if ((*prop_levelset2)[eid] > 0)
                    {
                        // index of the other fracture
                        pos = BaseLib::findIndex(br_matids, *frac2_matid);
                    }
                    assert(pos != std::numeric_limits<std::size_t>::max());
                    new_node_id = dup_newNodeIDs[pos];
                }
                else
                {
                    // junction nodes
                    const auto& jct_matids = getMaterialIdsForNode(
                        vec_junction_nodeID_matIDs, node_id);
                    auto const frac2_matid = BaseLib::findFirstNotEqualElement(
                        jct_matids, frac_matid);
                    assert(frac2_matid);
                    auto const frac2_id =
                        BaseLib::findIndex(vec_fracture_mat_IDs, *frac2_matid);
                    assert(frac2_id != std::numeric_limits<std::size_t>::max());
                    auto prop_levelset2 =
                        org_mesh.getProperties().getPropertyVector<double>(
                            "levelset" + std::to_string(frac2_id + 1));

                    //
                    if ((*prop_levelset2)[eid] <= 0)
                    {
                        // index of this frac
                        auto const pos =
                            BaseLib::findIndex(jct_matids, frac_matid);
                        assert(pos != std::numeric_limits<std::size_t>::max());
                        new_node_id = dup_newNodeIDs[pos];
                    }
                    else if ((*prop_levelset2)[eid] > 0)
                    {
                        // set the last duplicated node
                        new_node_id = dup_newNodeIDs.back();
                    }
                }

                // replace node
                e->setNode(i, new_nodes[new_node_id]);
            }
        }
    }

    // new mesh
    _output_mesh = std::make_unique<MeshLib::Mesh>(org_mesh.getName(),
                                                   new_nodes, new_eles);

    for (auto [name, property] : _org_mesh.getProperties())
    {
        if (auto p = dynamic_cast<MeshLib::PropertyVector<double>*>(property))
        {
            copyPropertyValues(*p, createProperty(*p));
        }
        else if (auto p =
                     dynamic_cast<MeshLib::PropertyVector<float>*>(property))
        {
            copyPropertyValues(*p, createProperty(*p));
        }
        else if (auto p = dynamic_cast<MeshLib::PropertyVector<int>*>(property))
        {
            copyPropertyValues(*p, createProperty(*p));
        }
        else if (auto p =
                     dynamic_cast<MeshLib::PropertyVector<unsigned>*>(property))
        {
            copyPropertyValues(*p, createProperty(*p));
        }
        else if (auto p =
                     dynamic_cast<MeshLib::PropertyVector<long>*>(property))
        {
            copyPropertyValues(*p, createProperty(*p));
        }
        else if (auto p = dynamic_cast<MeshLib::PropertyVector<long long>*>(
                     property))
        {
            copyPropertyValues(*p, createProperty(*p));
        }
        else if (auto p = dynamic_cast<MeshLib::PropertyVector<unsigned long>*>(
                     property))
        {
            copyPropertyValues(*p, createProperty(*p));
        }
        else if (auto p =
                     dynamic_cast<MeshLib::PropertyVector<unsigned long long>*>(
                         property))
        {
            copyPropertyValues(*p, createProperty(*p));
        }
        else if (auto p = dynamic_cast<MeshLib::PropertyVector<std::size_t>*>(
                     property))
        {
            copyPropertyValues(*p, createProperty(*p));
        }
        else if (auto p =
                     dynamic_cast<MeshLib::PropertyVector<char>*>(property))
        {
            copyPropertyValues(*p, createProperty(*p));
        }
        else
        {
            OGS_FATAL(
                "Mesh property '{:s}' of unhandled data type '{:s}'. Please "
                "check the data type of the mesh properties. The available "
                "data types are:"
                "\n\t double,"
                "\n\t float,"
                "\n\t int,"
                "\n\t unsigned,"
                "\n\t long,"
                "\n\t long long,"
                "\n\t unsigned long,"
                "\n\t unsigned long long,"
                "\n\t char.",
                property->getPropertyName(),
                typeid(*property).name());
        }
    }
    calculateTotalDisplacement(vec_vec_fracture_nodes.size(),
                               vec_junction_nodeID_matIDs.size());
}

template <typename T>
MeshLib::PropertyVector<T>* PostProcessTool::createProperty(
    MeshLib::PropertyVector<T> const& property)
{
    auto const item_type = property.getMeshItemType();
    auto const n_src_comp = property.getNumberOfGlobalComponents();
    // convert 2D vector to 3D. Otherwise Paraview Calculator filter does
    // not recognize it as a vector
    auto const n_dest_comp = (n_src_comp == 2) ? 3 : n_src_comp;

    auto new_property = MeshLib::getOrCreateMeshProperty<T>(
        *_output_mesh, property.getPropertyName(), item_type, n_dest_comp);

    if (item_type == MeshLib::MeshItemType::Node)
    {
        assert(new_property->size() ==
               _output_mesh->getNumberOfNodes() * n_dest_comp);
    }
    else if (item_type == MeshLib::MeshItemType::Cell)
    {
        assert(new_property->size() ==
               _output_mesh->getNumberOfElements() * n_dest_comp);
    }
    else
    {
        WARN(
            "Property '{:s}' cannot be created because its mesh item type "
            "'{:s}' is not supported.",
            property.getPropertyName(), toString(item_type));
        _output_mesh->getProperties().removePropertyVector(
            new_property->getPropertyName());
        return nullptr;
    }
    return new_property;
}

template <typename T>
void PostProcessTool::copyPropertyValues(
    MeshLib::PropertyVector<T> const& source_property,
    MeshLib::PropertyVector<T>* const destination_property)
{
    if (destination_property == nullptr)
    {
        // skip the copy, because the destination wasn't created.
        return;
    }

    auto const item_type = source_property.getMeshItemType();
    if (item_type == MeshLib::MeshItemType::Node)
    {
        auto const n_src_comp = source_property.getNumberOfGlobalComponents();
        auto const n_dest_comp =
            destination_property->getNumberOfGlobalComponents();
        // copy existing
        for (unsigned i = 0; i < _org_mesh.getNumberOfNodes(); i++)
        {
            auto last = std::copy_n(&source_property[i * n_src_comp],
                                    n_src_comp,
                                    &(*destination_property)[i * n_dest_comp]);
            // set zero for components not existing in the original
            std::fill_n(last, n_dest_comp - n_src_comp, T{0});
        }
        // copy duplicated
        for (auto itr : _map_dup_newNodeIDs)
        {
            for (unsigned k = 0; k < itr.second.size(); k++)
            {
                std::copy_n(
                    &(*destination_property)[itr.first * n_dest_comp],
                    n_dest_comp,
                    &(*destination_property)[itr.second[k] * n_dest_comp]);
            }
        }
    }
    else if (item_type == MeshLib::MeshItemType::Cell)
    {
        assert(source_property.size() == destination_property->size());
        std::copy(source_property.begin(), source_property.end(),
                  destination_property->begin());
    }
    else
    {
        OGS_FATAL(
            "Property '{:s}' values cannot be copied because its mesh item "
            "type '{:s}' is not supported. Unexpected error, because the "
            "destination property was created.",
            source_property.getPropertyName(), toString(item_type));
    }
}

void PostProcessTool::calculateTotalDisplacement(unsigned const n_fractures,
                                                 unsigned const n_junctions)
{
    auto const& u = *_output_mesh->getProperties().getPropertyVector<double>(
        "displacement");
    auto const n_u_comp = u.getNumberOfGlobalComponents();
    assert(u.size() == _output_mesh->getNodes().size() * 3);
    auto& total_u =
        *_output_mesh->getProperties().createNewPropertyVector<double>(
            "u", MeshLib::MeshItemType::Node, n_u_comp);
    total_u.resize(u.size());
    std::copy(cbegin(u), cend(u), begin(total_u));

    for (unsigned enrich_id = 0; enrich_id < n_fractures + n_junctions;
         enrich_id++)
    {
        // nodal value of levelset
        std::vector<double> nodal_levelset(
            _output_mesh->getNodes().size(),
            std::numeric_limits<double>::quiet_NaN());

        auto const& ele_levelset =
            *_output_mesh->getProperties().getPropertyVector<double>(
                "levelset" + std::to_string(enrich_id + 1));
        for (MeshLib::Element const* e : _output_mesh->getElements())
        {
            if (e->getDimension() != _output_mesh->getDimension())
            {
                continue;
            }
            const double e_levelset = ele_levelset[e->getID()];

            for (unsigned i = 0; i < e->getNumberOfNodes(); i++)
            {
                nodal_levelset[getNodeIndex(*e, i)] = e_levelset;
            }
        }

        // update total displacements
        auto const& g =
            *_output_mesh->getProperties().getPropertyVector<double>(
                "displacement_jump" + std::to_string(enrich_id + 1));
        for (unsigned i = 0; i < _output_mesh->getNodes().size(); i++)
        {
            for (int j = 0; j < n_u_comp; j++)
            {
                total_u[i * n_u_comp + j] +=
                    nodal_levelset[i] * g[i * n_u_comp + j];
            }
        }
    }
}

}  // namespace LIE
}  // namespace ProcessLib
