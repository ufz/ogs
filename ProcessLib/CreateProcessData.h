/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
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
        nonlinear_solvers);

}  // namespace ProcessLib
