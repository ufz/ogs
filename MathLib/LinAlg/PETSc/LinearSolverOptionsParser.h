/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MathLib/LinAlg/LinearSolverOptionsParser.h"
#include "MathLib/LinAlg/PETSc/PETScLinearSolver.h"

namespace BaseLib
{
class ConfigTree;
}

namespace MathLib
{
template <>
struct LinearSolverOptionsParser<PETScLinearSolver> final
{
    /// The method parses the linear solver options for PETSc.
    /// @param solver_prefix the prefix for the linear solver to distinguish
    /// different linear solvers for instance in the staggered schema
    /// @param config the part of the property tree (usually created from the
    /// linear solver section in the project file)
    /// @return the first item of the returned tuple is the solver prefix as
    /// string, the second item are all the options passed as string to PETSc
    /// options database
    std::tuple<std::string, std::string> parseNameAndOptions(
        std::string solver_prefix,
        BaseLib::ConfigTree const* const config) const;
};

}  // namespace MathLib
