/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHLIB_ODE_ODESOLVERBUILDER_H
#define MATHLIB_ODE_ODESOLVERBUILDER_H

#include <logog/include/logog.hpp>

#include "ODESolver.h"
#include "ConcreteODESolver.h"

#ifdef CVODE_FOUND
#include "CVodeSolver.h"
#endif

namespace BaseLib { class ConfigTree; }

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

    ERR("No ODE solver could be created. Maybe it is because you did not build"
        " OGS6 with support for any external ODE solver library.");
    std::abort();
}

//! @}

} // namespace ODE
} // namespace MathLib

#endif // MATHLIB_ODE_ODESOLVERBUILDER_H
