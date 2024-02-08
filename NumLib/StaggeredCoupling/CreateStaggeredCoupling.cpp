/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on November 21, 2023, 3:37 PM
 */

#include "CreateStaggeredCoupling.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "NumLib/ODESolver/ConvergenceCriterion.h"
#include "StaggeredCoupling.h"

namespace NumLib
{

/// The function returns a set of process names and the maximum iteration number
/// of a sub-coupling or nothing.
std::vector<LocalCouplingParameters> parseLocalCoupling(
    BaseLib::ConfigTree const& config, const std::size_t max_process_number)
{
    auto const local_coupling_configs =
        //! \ogs_file_param{prj__time_loop__global_process_coupling__local_coupling_processes}
        config.getConfigSubtreeList("local_coupling_processes");

    if (local_coupling_configs.empty())
    {
        return {};
    }

    std::vector<LocalCouplingParameters> all_local_coupling_parameters;
    std::vector<std::string> all_local_process_names;
    for (auto const& local_coupling_config : local_coupling_configs)
    {
        std::vector<std::string> process_names;

        for (
            auto name :
            local_coupling_config
                //! \ogs_file_param{prj__time_loop__global_process_coupling__local_coupling_processes__process_name}
                .getConfigParameterList<std::string>("process_name"))
        {
            if (std::find(process_names.begin(), process_names.end(), name) !=
                process_names.end())
            {
                OGS_FATAL(
                    "The name of locally coupled process, {}, is not unique.",
                    name);
            }
            process_names.push_back(name);

            all_local_process_names.push_back(name);
        }

        if (process_names.size() > max_process_number)
        {
            OGS_FATAL(
                "The number of the locally coupled processes is greater "
                "than the number of total coupled processes. "
                "Please check the number of elements in the tag "
                "'time_loop/global_process_coupling/"
                "local_coupling_processes' in the project file.");
        }

        INFO("There are {:d} locally coupled processes.", process_names.size());

        int max_iterations =
            local_coupling_config
                //! \ogs_file_param{prj__time_loop__global_process_coupling__local_coupling_processes__max_iter}
                .getConfigParameter<int>("max_iter");

        all_local_coupling_parameters.push_back(
            {process_names, max_iterations});
    }

    // std::adjacent_find only finds equal elements directly next to each other.
    // Therefore, a copy of the vector is sorted first and then it is checked
    // for duplicated element.
    std::vector<std::string> copy_all_local_process_names =
        all_local_process_names;
    std::sort(copy_all_local_process_names.begin(),
              copy_all_local_process_names.end());
    if (auto it = std::adjacent_find(copy_all_local_process_names.begin(),
                                     copy_all_local_process_names.end());
        it != copy_all_local_process_names.end())
    {
        OGS_FATAL(
            "There are process names appearing in multiple tags of "
            "'time_loop/global_process_coupling/local_coupling_processes'. For "
            "example, name {}",
            *it);
    }

    return all_local_coupling_parameters;
}

std::tuple<std::vector<std::unique_ptr<NumLib::ConvergenceCriterion>>,
           std::vector<LocalCouplingParameters>,
           int>
parseCoupling(BaseLib::ConfigTree const& config)
{
    auto const& coupling_config
        //! \ogs_file_param{prj__time_loop__global_process_coupling}
        = config.getConfigSubtreeOptional("global_process_coupling");

    std::vector<std::unique_ptr<NumLib::ConvergenceCriterion>>
        global_coupling_conv_criteria;

    std::vector<LocalCouplingParameters> all_local_coupling_parameters;

    int max_coupling_iterations = 1;
    if (coupling_config)
    {
        max_coupling_iterations
            //! \ogs_file_param{prj__time_loop__global_process_coupling__max_iter}
            = coupling_config->getConfigParameter<int>("max_iter");

        auto const& coupling_convergence_criteria_config =
            //! \ogs_file_param{prj__time_loop__global_process_coupling__convergence_criteria}
            coupling_config->getConfigSubtree("convergence_criteria");

        auto coupling_convergence_criterion_config =
            //! \ogs_file_param{prj__time_loop__global_process_coupling__convergence_criteria__convergence_criterion}
            coupling_convergence_criteria_config.getConfigSubtreeList(
                "convergence_criterion");
        std::transform(coupling_convergence_criterion_config.begin(),
                       coupling_convergence_criterion_config.end(),
                       std::back_inserter(global_coupling_conv_criteria),
                       [](BaseLib::ConfigTree const& c)
                       { return NumLib::createConvergenceCriterion(c); });

        all_local_coupling_parameters = parseLocalCoupling(
            *coupling_config, global_coupling_conv_criteria.size());
    }

    return {std::move(global_coupling_conv_criteria),
            std::move(all_local_coupling_parameters), max_coupling_iterations};
}

}  // namespace NumLib
