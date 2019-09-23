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

namespace NumLib
{

//! \addtogroup ODESolver
//! @{

//! Tag used to specify which nonlinear solver will be used.
enum class NonlinearSolverTag : bool {
    Picard /*!< Picard fixpoint iteration scheme */,
    Newton /*!< Newton-Raphson iteration scheme */
};

//! Tag used to specify the type of ODE.
enum class ODESystemTag : char
{
    /*! First order implicit quasi-linear ODE
     *
     * This is an ODE of the form
     * \f$ M(x,t)\cdot \dot x + K(x,t) \cdot x - b(x,t)
     *  =: r(\dot x, x, t) \stackrel{!}{=} 0 \f$
     */
    FirstOrderImplicitQuasilinear,
    NeumannBC // Sure, that's misuse of this enum, so sue me!
};

//! @}

} // namespace NumLib
