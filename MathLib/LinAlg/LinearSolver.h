/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHLIB_LINEAR_SOLVER_H
#define MATHLIB_LINEAR_SOLVER_H

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MathLib
{

// TODO also pass a name argument, or done within config?
/*! Creates a new linear solver instance.
 *
 * \tparam Solver the type of the linear solver to be created.
 *
 * \param config configuration options.
 */
template<typename Solver>
std::unique_ptr<Solver>
createLinearSolver(BaseLib::ConfigTree const*const config);

} // namespace MathLib

#endif // MATHLIB_LINEAR_SOLVER_H
