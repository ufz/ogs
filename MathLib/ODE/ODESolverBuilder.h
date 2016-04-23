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

#include "ODESolver.h"
#include "ConcreteODESolver.h"

#ifdef CVODE_FOUND
#include "CVodeSolver.h"
#endif

namespace BaseLib { class ConfigTree; }

namespace MathLib
{

template <unsigned NumEquations>
std::unique_ptr<ODESolver<NumEquations>> createODESolver(
    BaseLib::ConfigTree const& config)
{
#ifdef CVODE_FOUND
    return std::unique_ptr<ODESolver<NumEquations>>(
        new ConcreteODESolver<CVodeSolver, NumEquations>(config));
#else
    return nullptr;
#endif  // CVODE_FOUND
}

} // namespace MathLib

#endif // MATHLIB_ODE_ODESOLVERBUILDER_H
