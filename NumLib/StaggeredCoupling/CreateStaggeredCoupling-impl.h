/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on November 21, 2023, 5:12 PM
 */

#pragma once

#include <numeric>

#include "BaseLib/Error.h"
#include "CreateStaggeredCoupling.h"
#include "StaggeredCoupling.h"

namespace NumLib
{
using CouplingNodeVariant = std::variant<CouplingNode, RootCouplingNode>;

template <typename ProcessData>
void checkLocalCouplingParameters(
    std::vector<std::unique_ptr<ProcessData>> const& per_process_data,
    std::vector<LocalCouplingParameters> const& all_local_coupling_parameters)
{
    // Check whether the process names in the sub-coupling definitions exist in
    // the process data.
    for (auto const& local_coupling_parameters : all_local_coupling_parameters)
    {
        for (auto const& process_name : local_coupling_parameters.process_names)
        {
            if (std::none_of(
                    per_process_data.begin(),
                    per_process_data.end(),
                    [&process_name](auto const& process_data)
                    { return process_data->process_name == process_name; }))
            {
                OGS_FATAL(
                    "The given process name '{}' for the element "
                    "'time_loop/global_process_coupling/"
                    "local_coupling_processes/process_name' is not found "
                    "in the element 'time_loop/global_process_coupling/"
                    "local_coupling_processes/process_name' in the project "
                    "file.",
                    process_name);
            }
        }
    }
}

/// Create coupling nodes that do not have local-coupling nodes.
template <typename ProcessData>
std::vector<CouplingNodeVariant> createRegularCouplingNodes(
    std::vector<std::unique_ptr<ProcessData>> const& per_process_data,
    std::vector<LocalCouplingParameters> const& all_local_coupling_parameters,
    int const global_max_coupling_iterations,
    std::vector<std::unique_ptr<NumLib::ConvergenceCriterion>>&
        global_coupling_conv_criteria)
{
    std::vector<CouplingNodeVariant> coupling_nodes;

    for (auto const& process_data : per_process_data)
    {
        auto const& process_name = process_data->process_name;

        // If process_data->process_name occurs in local_coupling_parameters
        // do nothing
        if (std::any_of(
                all_local_coupling_parameters.begin(),
                all_local_coupling_parameters.end(),
                [&process_name](auto const& local_coupling_parameters)
                {
                    auto const& process_names =
                        local_coupling_parameters.process_names;
                    return std::find(process_names.begin(), process_names.end(),
                                     process_name) != process_names.end();
                }))
        {
            continue;
        }

        std::string const used_process_name =
            process_name.empty() ? "not given" : process_name;
        CouplingNode regular_node{
            used_process_name,
            std::move(global_coupling_conv_criteria[process_data->process_id]),
            global_max_coupling_iterations,
            process_data->process_id,
        };
        coupling_nodes.emplace_back(std::move(regular_node));
    }

    return coupling_nodes;
}

/// Create coupling nodes that have local-coupling nodes.
template <typename ProcessData>
std::vector<CouplingNodeVariant> createRootCouplingNodes(
    std::vector<std::unique_ptr<ProcessData>> const& per_process_data,
    std::vector<LocalCouplingParameters> const& all_local_coupling_parameters,
    int const global_max_coupling_iterations,
    std::vector<std::unique_ptr<NumLib::ConvergenceCriterion>>&
        global_coupling_conv_criteria)
{
    std::vector<CouplingNodeVariant> coupling_nodes;

    for (auto const& local_coupling_parameters : all_local_coupling_parameters)
    {
        RootCouplingNode root_node = {global_max_coupling_iterations, {}};

        for (auto const& local_process_name :
             local_coupling_parameters.process_names)
        {
            if (auto it = std::find_if(
                    per_process_data.begin(),
                    per_process_data.end(),
                    [&local_process_name](auto const& process_data) {
                        return process_data->process_name == local_process_name;
                    });
                it != per_process_data.end())
            {
                auto const& process_data = *it;

                CouplingNode regular_node{
                    process_data->process_name,
                    std::move(global_coupling_conv_criteria[process_data
                                                                ->process_id]),
                    local_coupling_parameters.max_iterations,
                    process_data->process_id};

                root_node.sub_coupling_nodes.emplace_back(
                    std::move(regular_node));
            }
        }
        coupling_nodes.emplace_back(std::move(root_node));
    }

    return coupling_nodes;
}

template <typename ProcessData>
std::vector<CouplingNodeVariant> createCouplingNodes(
    std::vector<std::unique_ptr<ProcessData>> const& per_process_data,
    std::vector<LocalCouplingParameters> const& all_local_coupling_parameters,
    int const global_max_coupling_iterations,
    std::vector<std::unique_ptr<NumLib::ConvergenceCriterion>>&
        global_coupling_conv_criteria)
{
    checkLocalCouplingParameters(per_process_data,
                                 all_local_coupling_parameters);

    // First, get the coupling nodes that do not have local-coupling nodes.
    std::vector<CouplingNodeVariant> coupling_nodes =
        createRegularCouplingNodes(
            per_process_data, all_local_coupling_parameters,
            global_max_coupling_iterations, global_coupling_conv_criteria);

    // Second, get the coupling nodes that have local-coupling nodes.
    std::vector<CouplingNodeVariant> root_coupling_nodes =
        createRootCouplingNodes(per_process_data, all_local_coupling_parameters,
                                global_max_coupling_iterations,
                                global_coupling_conv_criteria);

    std::size_t const num_coupling_nodes =
        coupling_nodes.size() +
        std::accumulate(
            root_coupling_nodes.begin(),
            root_coupling_nodes.end(),
            0,
            [](std::size_t accumulated_sizes, const auto& coupling_node)
            {
                return accumulated_sizes +
                       std::get<RootCouplingNode>(coupling_node)
                           .sub_coupling_nodes.size();
            });

    if (num_coupling_nodes != per_process_data.size())
    {
        OGS_FATAL(
            "The number of all coupling nodes including sub-nodes is not "
            "identical to the number of the processes! Please check the "
            "element by tag global_process_coupling in the project file.");
    }

    if (coupling_nodes.empty())
    {
        coupling_nodes = std::move(root_coupling_nodes);
    }
    else
    {
        coupling_nodes.reserve(coupling_nodes.size() +
                               root_coupling_nodes.size());

        std::move(std::begin(root_coupling_nodes),
                  std::end(root_coupling_nodes),
                  std::back_inserter(coupling_nodes));
    }

    return coupling_nodes;
}

/// Create a StaggeredCoupling instance from the given configuration.
template <typename ProcessData>
std::unique_ptr<StaggeredCoupling> createStaggeredCoupling(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ProcessData>> const& per_process_data)
{
    auto [global_coupling_conv_criteria, all_local_coupling_parameters,
          max_coupling_iterations] = parseCoupling(config);

    if (per_process_data.size() != global_coupling_conv_criteria.size())
    {
        OGS_FATAL(
            "The number of convergence criteria of the global "
            "staggered coupling loop is not identical to the number of the "
            "processes! Please check the element by tag "
            "global_process_coupling in the project file.");
    }

    auto coupling_nodes = createCouplingNodes(
        per_process_data, all_local_coupling_parameters,
        max_coupling_iterations, global_coupling_conv_criteria);

    return std::make_unique<StaggeredCoupling>(max_coupling_iterations,
                                               std::move(coupling_nodes));
}
}  // namespace NumLib
