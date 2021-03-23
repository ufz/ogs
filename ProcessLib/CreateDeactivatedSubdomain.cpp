/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#include "CreateDeactivatedSubdomain.h"

#include <Eigen/Dense>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "DeactivatedSubdomain.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
#include "MeshLib/Node.h"

namespace ProcessLib
{
template <typename IsActive>
static std::vector<MeshLib::Node*> extractInnerNodes(
    MeshLib::Mesh const& mesh,
    MeshLib::Mesh const& sub_mesh,
    IsActive is_active)
{
    auto* const bulk_node_ids =
        sub_mesh.getProperties().template getPropertyVector<std::size_t>(
            "bulk_node_ids", MeshLib::MeshItemType::Node, 1);
    if (bulk_node_ids == nullptr)
    {
        OGS_FATAL(
            "Bulk node ids map is not available in the deactivate subdomain "
            "mesh.");
    }

    std::vector<MeshLib::Node*> inner_nodes;
    // Reserve for more than needed, but not much, because almost all nodes will
    // be found and copied.
    inner_nodes.reserve(sub_mesh.getNumberOfNodes());
    std::copy_if(begin(sub_mesh.getNodes()), end(sub_mesh.getNodes()),
                 back_inserter(inner_nodes), [&](MeshLib::Node* const n) {
                     auto const bulk_node =
                         mesh.getNode((*bulk_node_ids)[n->getID()]);
                     const auto& connected_elements = bulk_node->getElements();

                     // Check whether this node is connected to an active
                     // element.
                     return std::all_of(begin(connected_elements),
                                        end(connected_elements), is_active);
                 });

    return inner_nodes;
}

static std::unique_ptr<DeactivatedSubdomainMesh> createDeactivatedSubdomainMesh(
    MeshLib::Mesh const& mesh, int const material_id)
{
    // An element is active if its material id matches the selected material id.
    auto is_active = [material_id, material_ids = *materialIDs(mesh)](
                         MeshLib::Element const* const e) {
        return material_id == material_ids[e->getID()];
    };

    auto const& elements = mesh.getElements();
    std::vector<MeshLib::Element*> deactivated_elements;
    std::copy_if(begin(elements), end(elements),
                 back_inserter(deactivated_elements),
                 [&](auto const e) { return is_active(e); });

    // Subdomain mesh consisting of deactivated elements.
    auto sub_mesh = MeshLib::createMeshFromElementSelection(
        "deactivate_subdomain_" + std::to_string(material_id),
        MeshLib::cloneElements(deactivated_elements));

    auto inner_nodes = extractInnerNodes(mesh, *sub_mesh, is_active);
    return std::make_unique<DeactivatedSubdomainMesh>(std::move(sub_mesh),
                                                      std::move(inner_nodes));
}

static MathLib::PiecewiseLinearInterpolation parseTimeIntervalOrCurve(
    std::optional<BaseLib::ConfigTree> const& time_interval_config,
    std::optional<std::string> const& curve_name,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    // Check for the error case first: only one of the configs must be used.
    if (time_interval_config && curve_name)
    {
        OGS_FATAL(
            "In the deactivate subdomain either a time interval or a curve "
            "must be given, not both.");
    }

    // Parse time interval and return a curve
    if (time_interval_config)
    {
        DBUG("Constructing time interval");
        auto const start_time =
            //! \ogs_file_param{prj__process_variables__process_variable__deactivated_subdomains__deactivated_subdomain__time_interval__start}
            time_interval_config->getConfigParameter<double>("start");

        auto const end_time =
            //! \ogs_file_param{prj__process_variables__process_variable__deactivated_subdomains__deactivated_subdomain__time_interval__end}
            time_interval_config->getConfigParameter<double>("end");

        // Using very large value for the curve's value, s.t. for any time from
        // start to end the whole subdomain is deactivated at once.
        return {{start_time, end_time},
                {std::numeric_limits<double>::max(),
                 std::numeric_limits<double>::max()},
                false};
    }

    // Try to find the curve.
    if (curve_name)
    {
        DBUG("Using curve '{:s}'", *curve_name);
        return *BaseLib::getOrError(curves, *curve_name,
                                    "Could not find curve.");
    }

    // If we got so far, there is an error: one of the configs must be
    // available.
    OGS_FATAL(
        "In the deactivate subdomain neither a time interval nor a curve are "
        "given. One of them must be specified.");
}

/// Returns a line segment represented by its begin and end points.
static std::pair<Eigen::Vector3d, Eigen::Vector3d> parseLineSegment(
    BaseLib::ConfigTree const& config)
{
    DBUG("Constructing line segment");
    auto const start =
        //! \ogs_file_param{prj__process_variables__process_variable__deactivated_subdomains__deactivated_subdomain__line_segment__start}
        config.getConfigParameter<std::vector<double>>("start");
    if (start.size() != 3)
    {
        OGS_FATAL(
            "For construction of a line segment the start point must be a "
            "3D point. Got a vector of size {}.",
            start.size());
    }

    auto const end =
        //! \ogs_file_param{prj__process_variables__process_variable__deactivated_subdomains__deactivated_subdomain__line_segment__end}
        config.getConfigParameter<std::vector<double>>("end");

    if (end.size() != 3)
    {
        OGS_FATAL(
            "For construction of a line segment the end point must be a "
            "3D point. Got a vector of size {}.",
            end.size());
    }
    return {Eigen::Vector3d{start[0], start[1], start[2]},
            Eigen::Vector3d{end[0], end[1], end[2]}};
}

std::unique_ptr<DeactivatedSubdomain const> createDeactivatedSubdomain(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& mesh,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    auto const& time_interval_config =
        //! \ogs_file_param{prj__process_variables__process_variable__deactivated_subdomains__deactivated_subdomain__time_interval}
        config.getConfigSubtreeOptional("time_interval");

    auto const& curve_name =
        //! \ogs_file_param{prj__process_variables__process_variable__deactivated_subdomains__deactivated_subdomain__time_curve}
        config.getConfigParameterOptional<std::string>("time_curve");
    auto const time_interval =
        parseTimeIntervalOrCurve(time_interval_config, curve_name, curves);

    auto const line_segment =
        //! \ogs_file_param{prj__process_variables__process_variable__deactivated_subdomains__deactivated_subdomain__line_segment}
        parseLineSegment(config.getConfigSubtree("line_segment"));

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

    auto const* const material_ids = materialIDs(mesh);
    if (material_ids == nullptr)
    {
        OGS_FATAL(
            "The mesh doesn't contain materialIDs for subdomain deactivation. "
            "The program terminates now.");
    }

    std::vector<std::unique_ptr<DeactivatedSubdomainMesh>>
        deactivated_subdomain_meshes;
    deactivated_subdomain_meshes.reserve(
        deactivated_subdomain_material_ids.size());

    for (int const id : deactivated_subdomain_material_ids)
    {
        deactivated_subdomain_meshes.push_back(
            createDeactivatedSubdomainMesh(mesh, id));
    }

    return std::make_unique<DeactivatedSubdomain const>(
        std::move(time_interval),
        line_segment,
        std::move(deactivated_subdomain_material_ids),
        std::move(deactivated_subdomain_meshes));
}

std::vector<std::unique_ptr<DeactivatedSubdomain const>>
createDeactivatedSubdomains(
    BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
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
                createDeactivatedSubdomain(subdomain_config, mesh, curves));
        }
    }
    return deactivated_subdomains;
}

}  // namespace ProcessLib
