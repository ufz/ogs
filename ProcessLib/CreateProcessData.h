/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessData.h"

namespace ProcessLib
{
std::vector<std::unique_ptr<ProcessData>> createPerProcessData(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<Process>> const& processes,
    std::map<std::string, std::unique_ptr<NumLib::NonlinearSolverBase>> const&
        nonlinear_solvers,
    bool const compensate_non_equilibrium_initial_residuum,
    std::vector<double> const& fixed_times_for_output,
    std::map<std::string, int>& local_coupling_processes);

}  // namespace ProcessLib
