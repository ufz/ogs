/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
    auto itr2 = std::find_if(
        vec_nodes.begin(), vec_nodes.end(),
        [&](MeshLib::Node const* node) { return node->getID() == node_id; });
    return (itr2 != vec_nodes.end());
}

std::vector<int> const& getMaterialIdsForNode(
    std::vector<std::pair<std::size_t, std::vector<int>>> const&
        vec_nodeID_matIDs,
    std::size_t nodeID)
{
    auto itr = std::find_if(
        vec_nodeID_matIDs.begin(), vec_nodeID_matIDs.end(),
        [&](std::pair<std::size_t, std::vector<int>> const& entry) {
            return entry.first == nodeID;
        });
    assert(itr != vec_nodeID_matIDs.end());
    return itr->second;
};

int findFirstNotEqualElement(std::vector<int> const& fracmatids, int frac1matid)
{
    auto itr_mat2 =
        std::find_if_not(fracmatids.begin(), fracmatids.end(),
                         [&](int matid) { return matid == frac1matid; });
    assert(itr_mat2 != fracmatids.end());
    return *itr_mat2;
};

int matid2fracid(std::vector<int> const& vec, int matid)
{
    auto itr = std::find(vec.begin(), vec.end(), matid);
    assert(itr != vec.end());
    return itr - vec.begin();
};

unsigned getpos_in_ids(std::vector<int> const& ids, int id)
{
    auto itr = std::find(ids.begin(), ids.end(), id);
    assert(itr != ids.end());
    auto const pos = itr - ids.begin();
    return pos;
}
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
                new MeshLib::Node(org_node->getCoords(), new_nodes.size());
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
            new MeshLib::Node(org_node->getCoords(), new_nodes.size());
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
                continue;

            auto const eid = org_e->getID();
            // keep original if the element has levelset=0
            if ((*prop_levelset)[eid] == 0)
                continue;

            // replace fracture nodes with duplicated ones
            MeshLib::Element* e = new_eles[eid];
            for (unsigned i = 0; i < e->getNumberOfNodes(); i++)
            {
                const auto node_id = e->getNodeIndex(i);
                if (!includesNodeID(vec_fracture_nodes, node_id))
                    continue;

                // list of duplicated node IDs
                auto itr = _map_dup_newNodeIDs.find(node_id);
                if (itr == _map_dup_newNodeIDs.end())
                    continue;
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
                    auto frac2_matid =
                        findFirstNotEqualElement(br_matids, frac_matid);
                    auto const frac2_id =
                        matid2fracid(vec_fracture_mat_IDs, frac2_matid);
                    auto prop_levelset2 =
                        org_mesh.getProperties().getPropertyVector<double>(
                            "levelset" + std::to_string(frac2_id + 1));
                    unsigned pos = 0;
                    if ((*prop_levelset2)[eid] == 0)
                    {
                        // index of this frac
                        pos = getpos_in_ids(br_matids, frac_matid);
                    }
                    else if ((*prop_levelset2)[eid] == 1)
                    {
                        // index of the other frac
                        pos = getpos_in_ids(br_matids, frac2_matid);
                    }
                    new_node_id = dup_newNodeIDs[pos];
                }
                else
                {
                    // junction nodes
                    const auto& jct_matids = getMaterialIdsForNode(
                        vec_junction_nodeID_matIDs, node_id);
                    auto frac2_matid =
                        findFirstNotEqualElement(jct_matids, frac_matid);
                    auto const frac2_id =
                        matid2fracid(vec_fracture_mat_IDs, frac2_matid);
                    auto prop_levelset2 =
                        org_mesh.getProperties().getPropertyVector<double>(
                            "levelset" + std::to_string(frac2_id + 1));

                    //
                    if ((*prop_levelset2)[eid] == 0)
                    {
                        // index of this frac
                        auto const pos = getpos_in_ids(jct_matids, frac_matid);
                        new_node_id = dup_newNodeIDs[pos];
                    }
                    else if ((*prop_levelset2)[eid] == 1)
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
    createProperties<int>();
    createProperties<double>();
    copyProperties<int>();
    copyProperties<double>();
    calculateTotalDisplacement(vec_vec_fracture_nodes.size(),
                               vec_junction_nodeID_matIDs.size());
}

template <typename T>
void PostProcessTool::createProperties()
{
    MeshLib::Properties const& src_properties = _org_mesh.getProperties();
    for (auto name : src_properties.getPropertyVectorNames())
    {
        if (!src_properties.existsPropertyVector<T>(name))
            continue;
        auto const* src_prop = src_properties.getPropertyVector<T>(name);

        auto const n_src_comp = src_prop->getNumberOfComponents();
        // convert 2D vector to 3D. Otherwise Paraview Calculator filter does
        // not recognize it as a vector
        auto const n_dest_comp = (n_src_comp == 2) ? 3 : n_src_comp;

        auto new_prop = MeshLib::getOrCreateMeshProperty<T>(
            *_output_mesh, name, src_prop->getMeshItemType(), n_dest_comp);

        if (src_prop->getMeshItemType() == MeshLib::MeshItemType::Node)
        {
            assert(new_prop->size() ==
                   _output_mesh->getNumberOfNodes() * n_dest_comp);
            (void)(new_prop);  // to avoid compilation warning.
        }
        else if (src_prop->getMeshItemType() == MeshLib::MeshItemType::Cell)
        {
            assert(new_prop->size() ==
                   _output_mesh->getNumberOfElements() * n_dest_comp);
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
        if (!src_properties.existsPropertyVector<T>(name))
            continue;
        auto const* src_prop = src_properties.getPropertyVector<T>(name);
        auto* dest_prop =
            _output_mesh->getProperties().getPropertyVector<T>(name);

        if (src_prop->getMeshItemType() == MeshLib::MeshItemType::Node)
        {
            auto const n_src_comp = src_prop->getNumberOfComponents();
            auto const n_dest_comp = dest_prop->getNumberOfComponents();
            // copy existing
            for (unsigned i = 0; i < _org_mesh.getNumberOfNodes(); i++)
            {
                for (int j = 0; j < n_src_comp; j++)
                    (*dest_prop)[i * n_dest_comp + j] =
                        (*src_prop)[i * n_src_comp + j];
                // set zero for components not existing in the original
                for (int j = n_src_comp; j < n_dest_comp; j++)
                    (*dest_prop)[i * n_dest_comp + j] = 0;
            }
            // copy duplicated
            for (auto itr : _map_dup_newNodeIDs)
            {
                for (int j = 0; j < n_dest_comp; j++)
                    for (unsigned k = 0; k < itr.second.size(); k++)
                        (*dest_prop)[itr.second[k] * n_dest_comp + j] =
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

void PostProcessTool::calculateTotalDisplacement(unsigned const n_fractures,
                                                 unsigned const n_junctions)
{
    auto const& u = *_output_mesh->getProperties().getPropertyVector<double>(
        "displacement");
    auto const n_u_comp = u.getNumberOfComponents();
    assert(u.size() == _output_mesh->getNodes().size() * 3);
    auto& total_u =
        *_output_mesh->getProperties().createNewPropertyVector<double>(
            "u", MeshLib::MeshItemType::Node, n_u_comp);
    total_u.resize(u.size());
    for (unsigned i = 0; i < _output_mesh->getNodes().size(); i++)
    {
        for (int j = 0; j < n_u_comp; j++)
            total_u[i * n_u_comp + j] = u[i * n_u_comp + j];
    }

    for (unsigned enrich_id = 0; enrich_id < n_fractures + n_junctions;
         enrich_id++)
    {
        // nodal value of levelset
        std::vector<double> nodal_levelset(_output_mesh->getNodes().size(),
                                           0.0);
        auto const& ele_levelset =
            *_output_mesh->getProperties().getPropertyVector<double>(
                "levelset" + std::to_string(enrich_id + 1));
        for (MeshLib::Element const* e : _output_mesh->getElements())
        {
            if (e->getDimension() != _output_mesh->getDimension())
                continue;
            const double e_levelset = ele_levelset[e->getID()];
            if (e_levelset == 0)
                continue;

            for (unsigned i = 0; i < e->getNumberOfNodes(); i++)
                nodal_levelset[e->getNodeIndex(i)] = e_levelset;
        }

        // update total displacements
        auto const& g =
            *_output_mesh->getProperties().getPropertyVector<double>(
                "displacement_jump" + std::to_string(enrich_id + 1));
        for (unsigned i = 0; i < _output_mesh->getNodes().size(); i++)
        {
            for (int j = 0; j < n_u_comp; j++)
                total_u[i * n_u_comp + j] +=
                    nodal_levelset[i] * g[i * n_u_comp + j];
        }
    }
}

}  // namespace LIE
}  // namespace ProcessLib
