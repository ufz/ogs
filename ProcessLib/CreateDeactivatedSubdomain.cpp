/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#include "CreateDeactivatedSubdomain.h"

#include <Eigen/Core>
#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/algorithm/copy_if.hpp>
#include <range/v3/algorithm/partition_copy.hpp>
#include <range/v3/range/conversion.hpp>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "DeactivatedSubdomain.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Utils/DuplicateMeshComponents.h"
#include "MeshLib/Utils/createMeshFromElementSelection.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/Utils.h"

namespace ProcessLib
{
template <typename IsActive>
static std::pair<std::vector<std::size_t>, std::vector<std::size_t>>
extractInnerAndOuterNodes(MeshLib::Mesh const& mesh,
                          MeshLib::Mesh const& sub_mesh,
                          IsActive is_active)
{
    auto* const bulk_node_ids = MeshLib::bulkNodeIDs(sub_mesh);
    if (bulk_node_ids == nullptr)
    {
        OGS_FATAL(
            "Bulk node ids map is not available in the deactivate subdomain "
            "mesh.");
    }

    std::vector<std::size_t> inner_nodes;
    // Reserve for more than needed, but not much, because almost all nodes will
    // be the inner nodes.
    inner_nodes.reserve(sub_mesh.getNumberOfNodes());

    std::vector<std::size_t> outer_nodes;

    ranges::partition_copy(
        sub_mesh.getNodes() | MeshLib::views::ids,
        back_inserter(inner_nodes),
        back_inserter(outer_nodes),
        [&](std::size_t const n)
        {
            auto const bulk_node = mesh.getNode((*bulk_node_ids)[n]);
            const auto& connected_elements =
                mesh.getElementsConnectedToNode(*bulk_node);
            // Check whether this node is connected only to an active
            // elements. Then it is an inner node, and outer otherwise.
            return ranges::all_of(connected_elements | MeshLib::views::ids,
                                  is_active);
        });

    return {std::move(inner_nodes), std::move(outer_nodes)};
}

static std::vector<std::vector<std::size_t>> extractElementsAlongOuterNodes(
    MeshLib::Mesh const& mesh, MeshLib::Mesh const& sub_mesh,
    std::vector<std::size_t> const& outer_nodes)
{
    auto const& bulk_node_ids =
        *sub_mesh.getProperties().template getPropertyVector<std::size_t>(
            MeshLib::getBulkIDString(MeshLib::MeshItemType::Node),
            MeshLib::MeshItemType::Node, 1);

    auto to_bulk_node_id =
        ranges::views::transform([&bulk_node_ids](std::size_t const node_id)
                                 { return bulk_node_ids[node_id]; });
    auto connected_elements = ranges::views::transform(
        [&mesh](std::size_t const node_id)
        {
            return mesh.getElementsConnectedToNode(node_id) |
                   MeshLib::views::ids | ranges::to<std::vector>;
        });

    return outer_nodes | to_bulk_node_id | connected_elements |
           ranges::to<std::vector>;
}

static DeactivatedSubdomainMesh createDeactivatedSubdomainMesh(
    MeshLib::Mesh const& mesh,
    std::unordered_set<int> const& deactivated_subdomain_material_ids)
{
    // An element is active if its material id matches one of the deactivated
    // subdomain material ids.
    auto is_active =
        [&deactivated_subdomain_material_ids,
         &material_ids = *materialIDs(mesh)](std::size_t element_id)
    {
        return deactivated_subdomain_material_ids.contains(
            material_ids[element_id]);
    };

    auto const& elements = mesh.getElements();
    std::vector<MeshLib::Element*> deactivated_elements;
    ranges::copy_if(elements, back_inserter(deactivated_elements), is_active,
                    [](auto const e) { return e->getID(); });

    auto bulk_element_ids = deactivated_elements | MeshLib::views::ids |
                            ranges::to<std::unordered_set>;

    static int mesh_number = 0;
    // Subdomain mesh consisting of deactivated elements.
    auto sub_mesh = MeshLib::createMeshFromElementSelection(
        "deactivate_subdomain_" + std::to_string(mesh_number++),
        MeshLib::cloneElements(deactivated_elements));

    auto [inner_nodes, outer_nodes] =
        extractInnerAndOuterNodes(mesh, *sub_mesh, is_active);

    auto outer_nodes_elements =
        extractElementsAlongOuterNodes(mesh, *sub_mesh, outer_nodes);

    return {std::move(*sub_mesh), std::move(bulk_element_ids),
            std::move(inner_nodes), std::move(outer_nodes),
            std::move(outer_nodes_elements)};
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
        // Return a copy of the curve because the time interval in the other
        // branch returns a temporary.
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
            "For construction of a line segment the start point must be a 3D "
            "point. Got a vector of size {}.",
            start.size());
    }

    auto const end =
        //! \ogs_file_param{prj__process_variables__process_variable__deactivated_subdomains__deactivated_subdomain__line_segment__end}
        config.getConfigParameter<std::vector<double>>("end");

    if (end.size() != 3)
    {
        OGS_FATAL(
            "For construction of a line segment the end point must be a 3D "
            "point. Got a vector of size {}.",
            end.size());
    }
    return {Eigen::Vector3d{start[0], start[1], start[2]},
            Eigen::Vector3d{end[0], end[1], end[2]}};
}

DeactivatedSubdomain createDeactivatedSubdomain(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& mesh,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
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
    auto time_interval =
        parseTimeIntervalOrCurve(time_interval_config, curve_name, curves);

    auto const line_segment_config =
        //! \ogs_file_param{prj__process_variables__process_variable__deactivated_subdomains__deactivated_subdomain__line_segment}
        config.getConfigSubtreeOptional("line_segment");

    if (time_interval_config && line_segment_config)
    {
        OGS_FATAL(
            "When using time interval for subdomain deactivation a line "
            "segment must not be specified.");
    }

    if (curve_name && !line_segment_config)
    {
        OGS_FATAL(
            "When using curve for subdomain deactivation a line segment must "
            "be specified.");
    }

    // If time interval was specified then an empty optional line segment is
    // used *internally* because the whole selected material ids subdomain will
    // be deactivated.
    std::optional<std::pair<Eigen::Vector3d, Eigen::Vector3d>> line_segment;
    if (curve_name)
    {
        line_segment = parseLineSegment(*line_segment_config);
    }

    ParameterLib::Parameter<double>* boundary_value_parameter = nullptr;
    auto boundary_value_parameter_name =
        //! \ogs_file_param{prj__process_variables__process_variable__deactivated_subdomains__deactivated_subdomain__boundary_parameter}
        config.getConfigParameterOptional<std::string>("boundary_parameter");
    if (boundary_value_parameter_name)
    {
        DBUG("Using parameter {:s}", *boundary_value_parameter_name);
        boundary_value_parameter = &ParameterLib::findParameter<double>(
            *boundary_value_parameter_name, parameters, 1, &mesh);
    }

    auto deactivated_subdomain_material_ids =
        //! \ogs_file_param{prj__process_variables__process_variable__deactivated_subdomains__deactivated_subdomain__material_ids}
        config.getConfigParameter("material_ids", std::vector<int>{}) |
        ranges::to<std::unordered_set>();

    if (deactivated_subdomain_material_ids.empty())
    {
        OGS_FATAL(
            "The material IDs of the deactivated subdomains are not given. The "
            "program terminates now.");
    }

    if (materialIDs(mesh) == nullptr)
    {
        OGS_FATAL(
            "The mesh doesn't contain materialIDs for subdomain deactivation. "
            "The program terminates now.");
    }

    auto deactivated_subdomain_mesh = createDeactivatedSubdomainMesh(
        mesh, deactivated_subdomain_material_ids);

    return {std::move(time_interval), line_segment,
            std::move(deactivated_subdomain_mesh), boundary_value_parameter};
}

std::vector<DeactivatedSubdomain> createDeactivatedSubdomains(
    BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    std::vector<DeactivatedSubdomain> deactivated_subdomains;
    // Deactivated subdomains
    if (auto subdomains_config =
            //! \ogs_file_param{prj__process_variables__process_variable__deactivated_subdomains}
        config.getConfigSubtreeOptional("deactivated_subdomains"))
    {
        INFO("There are subdomains being deactivated.");

        auto const deactivated_subdomain_configs =
            //! \ogs_file_param{prj__process_variables__process_variable__deactivated_subdomains__deactivated_subdomain}
            subdomains_config->getConfigSubtreeList("deactivated_subdomain");
        std::transform(std::begin(deactivated_subdomain_configs),
                       std::end(deactivated_subdomain_configs),
                       std::back_inserter(deactivated_subdomains),
                       [&](auto const& config) {
                           return createDeactivatedSubdomain(
                               config, mesh, parameters, curves);
                       });
    }
    return deactivated_subdomains;
}

}  // namespace ProcessLib
