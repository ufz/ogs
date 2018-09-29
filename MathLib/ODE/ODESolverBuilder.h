/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <logog/include/logog.hpp>

#include "BaseLib/Error.h"
#include "ODESolver.h"
#include "ConcreteODESolver.h"

#ifdef CVODE_FOUND
#include "CVodeSolver.h"
#endif

namespace BaseLib
{
class ConfigTree;
}

namespace MathLib
{
namespace ODE
{
//! \addtogroup ExternalODESolverInterface
//! @{

/*! Creates a new ODESolver instance from the given \c config.
 *
 * \tparam NumEquations the number of equations in the ODE system to be solved.
 */
template <unsigned NumEquations>
std::unique_ptr<ODESolver<NumEquations>> createODESolver(
    BaseLib::ConfigTree const& config)
{
#ifdef CVODE_FOUND
    return std::unique_ptr<ODESolver<NumEquations>>(
        new ConcreteODESolver<CVodeSolver, NumEquations>(config));
#endif
    (void)config;  // Unused parameter warning if no library is available.

    OGS_FATAL(
        "No ODE solver could be created. Maybe it is because you did not build"
        " OGS6 with support for any external ODE solver library.");
}

//! @}

}  // namespace ODE
}  // namespace MathLib
