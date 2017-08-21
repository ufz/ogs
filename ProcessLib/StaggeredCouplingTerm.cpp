/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   StaggeredCouplingTerm.cpp
 *
 * Created on November 7, 2016, 12:14 PM
 */

#include "StaggeredCouplingTerm.h"

#include "MathLib/LinAlg/LinAlg.h"
#include "Process.h"

namespace ProcessLib
{

StaggeredCouplingTerm::StaggeredCouplingTerm(
    std::unordered_map<std::type_index, Process const&> const&
        coupled_processes_,
    std::unordered_map<std::type_index, GlobalVector const&> const& coupled_xs_,
    const double dt_, const bool empty_)
    : coupled_processes(coupled_processes_),
      coupled_xs(coupled_xs_),
      dt(dt_),
      empty(empty_)
{
    for (auto const& coupled_x_pair : coupled_xs)
    {
        auto const& coupled_x = coupled_x_pair.second;
        MathLib::LinAlg::setLocalAccessibleVector(coupled_x);
    }

    for (auto const& coupled_process_pair : coupled_processes)
    {
        auto const& coupled_pcs = coupled_process_pair.second;
        auto const prevous_time_x = coupled_pcs.getPreviousTimeStepSolution();
        if (prevous_time_x)
        {
            MathLib::LinAlg::setLocalAccessibleVector(*prevous_time_x);
        }
    }
}

std::unordered_map<std::type_index, const std::vector<double>>
getCurrentLocalSolutionsOfCoupledProcesses(
    const std::unordered_map<std::type_index, GlobalVector const&>&
        global_coupled_xs,
    const std::vector<GlobalIndexType>& indices)
{
    std::unordered_map<std::type_index, const std::vector<double>>
        local_coupled_xs;

    // Get local nodal solutions of the coupled equations.
    for (auto const& global_coupled_x_pair : global_coupled_xs)
    {
        auto const& coupled_x = global_coupled_x_pair.second;
        auto const local_coupled_x = coupled_x.get(indices);
        BaseLib::insertIfTypeIndexKeyUniqueElseError(
            local_coupled_xs, global_coupled_x_pair.first, local_coupled_x,
            "local_coupled_x");
    }
    return local_coupled_xs;
}

}  // end of ProcessLib
