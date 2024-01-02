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

namespace BaseLib
{
class ConfigTree;
}
#include "MathLib/LinAlg/EigenLis/EigenLisLinearSolver.h"
#include "MathLib/LinAlg/LinearSolverOptionsParser.h"

namespace MathLib
{
template <>
struct LinearSolverOptionsParser<EigenLisLinearSolver> final
{
    /// The method parses the linear solver options for LIS.
    /// @param prefix the prefix for the linear solver to distinguish different
    /// linear solvers for instance in the staggered schema
    /// @param config the part of the property tree (usually created from the
    /// linear solver section in the project file)
    /// @return the first item of the returned tuple is the solver prefix as
    /// string, the second item are all the options as a string that will be
    /// passed to LIS via lis_solver_set_options()
    std::tuple<std::string, std::string> parseNameAndOptions(
        std::string const& prefix,
        BaseLib::ConfigTree const* const config) const;
};

}  // namespace MathLib
