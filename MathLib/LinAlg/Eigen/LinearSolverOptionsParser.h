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

#include "EigenLinearSolver.h"
#include "EigenOption.h"
#include "MathLib/LinAlg/LinearSolverOptionsParser.h"

namespace BaseLib
{
class ConfigTree;
}

namespace MathLib
{
template <>
struct LinearSolverOptionsParser<EigenLinearSolver> final
{
    /// The method parses the linear solver options for Eigen.
    /// @param prefix the prefix for the linear solver to distinguish different
    /// linear solvers for instance in the staggered schema
    /// @param solver_config the part of the property tree (usually created from
    /// the linear solver section in the project file)
    /// @return the first item of the returned tuple is the solver prefix as
    /// string, the second item are all the options as an EigenOption object
    std::tuple<std::string, EigenOption> parseNameAndOptions(
        std::string const& prefix,
        BaseLib::ConfigTree const* const solver_config) const;
};

}  // namespace MathLib
