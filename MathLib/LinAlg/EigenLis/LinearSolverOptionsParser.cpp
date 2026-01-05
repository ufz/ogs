// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "LinearSolverOptionsParser.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Logging.h"
#include "MathLib/LinAlg/LinearSolverOptions.h"

namespace MathLib
{
std::tuple<std::string, std::string>
LinearSolverOptionsParser<EigenLisLinearSolver>::parseNameAndOptions(
    std::string const& prefix, BaseLib::ConfigTree const* const config) const
{
    std::string lis_options = "-initx_zeros 0";

    if (config)
    {
        ignoreOtherLinearSolvers(*config, "lis");
        //! \ogs_file_param{prj__linear_solvers__linear_solver__lis}
        if (auto s = config->getConfigParameterOptional<std::string>("lis"))
        {
            if (!s->empty())
            {
                lis_options += " " + *s;
                INFO("Lis options: '{:s}'", lis_options);
            }
        }
    }
    return {prefix, lis_options};
}
}  // namespace MathLib
