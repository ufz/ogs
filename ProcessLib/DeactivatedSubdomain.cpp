/**
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * File:   DeactivatedSubdomain.cpp
 *
 * Created on November 29, 2018, 10:50 AM
 */
#include "DeactivatedSubdomain.h"

#include <logog/include/logog.hpp>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"

namespace ProcessLib
{
const std::string DeactivatedSubdomain::name_of_paramater_of_zero =
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
        deactivated_sudomain_meshes_)
    : time_interval(std::move(time_interval_)),
      materialIDs(std::move(materialIDs_)),
      deactivated_sudomain_meshes(std::move(deactivated_sudomain_meshes_))
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

    std::vector<int> deactivated_subdomain_material_ids;
    deactivated_subdomain_material_ids =
        //! \ogs_file_param{prj__process_variables__process_variable__deactivated_subdomains__deactivated_subdomain__material_ids}
        config.getConfigParameter<std::vector<int>>("material_ids",
                                                    std::vector<int>{});
    if (!deactivated_subdomain_material_ids.empty())
    {
        std::sort(deactivated_subdomain_material_ids.begin(),
                  deactivated_subdomain_material_ids.end());

        auto const* const material_ids = MeshLib::materialIDs(mesh);

        std::vector<std::unique_ptr<DeactivetedSubdomainMesh>>
            deactivated_sudomain_meshes;
        deactivated_sudomain_meshes.reserve(
            deactivated_subdomain_material_ids.size());

        for (auto const ids : deactivated_subdomain_material_ids)
        {
            std::vector<MeshLib::Element*> deactivated_elements;
            std::vector<MeshLib::Node*> deactivated_nodes;

            // temporary vector to enhance node searching.
            std::vector<bool> deactivation_flag_of_nodes(
                mesh.getNumberOfNodes(), false);

            for (std::size_t i = 0; i < mesh.getNumberOfElements(); i++)
            {
                if (ids != (*material_ids)[i])
                    continue;

                auto* element = mesh.getElement(i);
                deactivated_elements.push_back(
                    const_cast<MeshLib::Element*>(element));

                for (unsigned i = 0; i < element->getNumberOfNodes(); i++)
                {
                    auto const* const node = element->getNode(i);
                    const auto& connected_elements = node->getElements();

                    if (deactivation_flag_of_nodes[node->getID()])
                        continue;

                    // Check whether this node is in an activated element.
                    if (std::find_if(
                            connected_elements.begin(),
                            connected_elements.end(),
                            [&](auto const* const connected_elem) -> bool {
                                return ids !=
                                       (*material_ids)[connected_elem->getID()];
                            }) != connected_elements.end())
                        continue;

                    deactivated_nodes.push_back(
                        const_cast<MeshLib::Node*>(node));
                    deactivation_flag_of_nodes[node->getID()] = true;
                }
            }

            auto bc_mesh = MeshLib::createMeshFromElementSelection(
                "deactivate_subdomain" + std::to_string(ids),
                MeshLib::cloneElements(deactivated_elements));

            deactivated_sudomain_meshes.emplace_back(
                std::make_unique<DeactivetedSubdomainMesh>(
                    std::move(bc_mesh), std::move(deactivated_nodes)));
        }

        return std::make_unique<DeactivatedSubdomain const>(
            std::move(time_interval),
            std::move(deactivated_subdomain_material_ids),
            std::move(deactivated_sudomain_meshes));
    }
    else
    {
        OGS_FATAL(
            "The material IDs of the deactivated subdomains are not "
            "given.");
    }

    return nullptr;
}

std::vector<std::unique_ptr<DeactivatedSubdomain const>>
createDeactivatedSubdomains(BaseLib::ConfigTree const& config,
                            MeshLib::Mesh const& mesh)
{
    std::vector<std::unique_ptr<DeactivatedSubdomain const>>
        deactivated_subdomains;
    // Deactivated subdomains
    //! \ogs_file_param{prj__process_variables__process_variable__deactivated_subdomains}
    if (auto subdomains_config =
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
