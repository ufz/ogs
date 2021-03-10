/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * File:   DeactivatedSubdomain.cpp
 *
 * Created on November 29, 2018, 10:50 AM
 */
#include "DeactivatedSubdomain.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
#include "MeshLib/Node.h"

namespace ProcessLib
{
const std::string DeactivatedSubdomain::zero_parameter_name =
    "zero_for_element_deactivation_approach";

DeactivetedSubdomainMesh::DeactivetedSubdomainMesh(
    std::unique_ptr<MeshLib::Mesh> deactivated_subdomain_mesh_,
    std::vector<MeshLib::Node*>&& inactive_nodes_)
    : mesh(std::move(deactivated_subdomain_mesh_)),
      inactive_nodes(std::move(inactive_nodes_))
{
}

DeactivatedSubdomain::DeactivatedSubdomain(
    std::unique_ptr<BaseLib::TimeInterval> time_interval_,
    std::vector<int>&& materialIDs_,
    std::vector<std::unique_ptr<DeactivetedSubdomainMesh>>&&
        deactivated_subdomain_meshes_)
    : time_interval(std::move(time_interval_)),
      materialIDs(std::move(materialIDs_)),
      deactivated_subdomain_meshes(std::move(deactivated_subdomain_meshes_))
{
}

bool DeactivatedSubdomain::includesTimeOf(double const t) const
{
    return time_interval->contains(t);
}

std::unique_ptr<DeactivatedSubdomain const> createDeactivatedSubdomain(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& mesh)
{
    //! \ogs_file_param{prj__process_variables__process_variable__deactivated_subdomains__deactivated_subdomain__time_interval}
    config.peekConfigParameter<std::string>("time_interval");
    auto time_interval = BaseLib::createTimeInterval(config);

    auto deactivated_subdomain_material_ids =
        //! \ogs_file_param{prj__process_variables__process_variable__deactivated_subdomains__deactivated_subdomain__material_ids}
        config.getConfigParameter("material_ids", std::vector<int>{});

    if (deactivated_subdomain_material_ids.empty())
    {
        OGS_FATAL(
            "The material IDs of the deactivated subdomains are not given. The "
            "program terminates now.");
    }

    std::sort(deactivated_subdomain_material_ids.begin(),
              deactivated_subdomain_material_ids.end());

    auto const* const material_ids = MeshLib::materialIDs(mesh);
    if (material_ids == nullptr)
    {
        OGS_FATAL(
            "The mesh doesn't contain materialIDs for subdomain deactivation. "
            "The program terminates now.");
    }

    std::vector<std::unique_ptr<DeactivetedSubdomainMesh>>
        deactivated_subdomain_meshes;
    deactivated_subdomain_meshes.reserve(
        deactivated_subdomain_material_ids.size());

    for (auto const ids : deactivated_subdomain_material_ids)
    {
        auto const& nodes = mesh.getNodes();
        std::vector<std::size_t> deactivated_bulk_node_ids;
        for (auto const& node : nodes)
        {
            const auto& connected_elements = node->getElements();

            // Check whether this node is in an activated element.
            if (std::find_if(
                    connected_elements.begin(),
                    connected_elements.end(),
                    [&](auto const* const connected_elem) -> bool {
                        return ids != (*material_ids)[connected_elem->getID()];
                    }) != connected_elements.end())
            {
                continue;
            }

            deactivated_bulk_node_ids.push_back(node->getID());
        }

        auto const& elements = mesh.getElements();
        std::vector<MeshLib::Element*> deactivated_elements;
        for (auto const& element : elements)
        {
            if (ids != (*material_ids)[element->getID()])
            {
                continue;
            }

            deactivated_elements.push_back(
                const_cast<MeshLib::Element*>(element));
        }

        auto bc_mesh = MeshLib::createMeshFromElementSelection(
            "deactivate_subdomain" + std::to_string(ids),
            MeshLib::cloneElements(deactivated_elements));

        auto const& new_mesh_properties = bc_mesh->getProperties();
        if (!new_mesh_properties.template existsPropertyVector<std::size_t>(
                "bulk_node_ids"))
        {
            OGS_FATAL(
                "Bulk node ids map expected in the construction of the mesh "
                "subset.");
        }
        auto const& bulk_node_ids_map =
            *new_mesh_properties.template getPropertyVector<std::size_t>(
                "bulk_node_ids", MeshLib::MeshItemType::Node, 1);

        std::vector<MeshLib::Node*> deactivated_nodes;
        auto const& nodes_in_bc_mesh = bc_mesh->getNodes();
        for (std::size_t i = 0; i < bulk_node_ids_map.size(); i++)
        {
            auto const found_iterator = std::find(
                deactivated_bulk_node_ids.begin(),
                deactivated_bulk_node_ids.end(), bulk_node_ids_map[i]);

            if (std::end(deactivated_bulk_node_ids) == found_iterator)
            {
                continue;
            }
            deactivated_nodes.push_back(
                const_cast<MeshLib::Node*>(nodes_in_bc_mesh[i]));

            deactivated_bulk_node_ids.erase(found_iterator);
        }

        deactivated_subdomain_meshes.emplace_back(
            std::make_unique<DeactivetedSubdomainMesh>(
                std::move(bc_mesh), std::move(deactivated_nodes)));
    }

    return std::make_unique<DeactivatedSubdomain const>(
        std::move(time_interval),
        std::move(deactivated_subdomain_material_ids),
        std::move(deactivated_subdomain_meshes));
}

std::vector<std::unique_ptr<DeactivatedSubdomain const>>
createDeactivatedSubdomains(BaseLib::ConfigTree const& config,
                            MeshLib::Mesh const& mesh)
{
    std::vector<std::unique_ptr<DeactivatedSubdomain const>>
        deactivated_subdomains;
    // Deactivated subdomains
    if (auto subdomains_config =
            //! \ogs_file_param{prj__process_variables__process_variable__deactivated_subdomains}
        config.getConfigSubtreeOptional("deactivated_subdomains"))
    {
        INFO("There are subdomains being deactivated.");

        for (
            auto subdomain_config :
            //! \ogs_file_param{prj__process_variables__process_variable__deactivated_subdomains__deactivated_subdomain}
            subdomains_config->getConfigSubtreeList("deactivated_subdomain"))
        {
            deactivated_subdomains.emplace_back(
                createDeactivatedSubdomain(subdomain_config, mesh));
        }
    }
    return deactivated_subdomains;
}

}  // namespace ProcessLib
