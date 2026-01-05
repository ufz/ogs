// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "LinearSolverOptions.h"

#include <set>

#include "BaseLib/ConfigTree.h"

//! Configuration tag names of all known linear solvers for their
//! configuration in the project file.
//! Add your tag name here when you add a new solver.
static std::set<std::string> known_linear_solvers{"eigen", "lis", "petsc"};

namespace MathLib
{
void ignoreOtherLinearSolvers(const BaseLib::ConfigTree& config,
                              const std::string& solver_name)
{
    for (auto const& s : known_linear_solvers)
    {
        if (s != solver_name)
        {
            config.ignoreConfigParameter(s);
        }
    }
}

}  // namespace MathLib
