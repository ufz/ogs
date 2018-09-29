/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/DOF/MatrixProviderUser.h"

namespace NumLib
{
//! \addtogroup ODESolver
//! @{

//! Status flags telling the NonlinearSolver if an iteration succeeded.
enum class IterationResult : char
{
    SUCCESS,
    FAILURE,
    REPEAT_ITERATION
};

/*! Collection of basic methods every equation system must provide.
 *
 */
class EquationSystem : public NumLib::MatrixSpecificationsProvider
{
public:
    /*! Check whether this is actually a linear equation system.
     *
     * \remark
     * Depending on its parameters an in general nonlinear equation system
     * can be linear in special cases. With this method it is possible to
     * detect that at runtime and thus save some computations.
     */
    virtual bool isLinear() const = 0;

    /*! Prepares a new iteration in the solution process of this equation.
     *
     * \param iter the current iteration number, starting from 1.
     * \param x    the current approximate solution of the equation.
     */
    virtual void preIteration(const unsigned iter, GlobalVector const& x)
    {
        (void)iter;
        (void)x;  // by default do nothing
    }

    /*! Post-processes an iteration in the solution process of this equation.
     *
     * \param x the current approximate solution of the equation.
     *
     * \return A status flag indicating id the current iteration succeeded.
     */
    virtual IterationResult postIteration(GlobalVector const& x)
    {
        (void)x;  // by default do nothing
        return IterationResult::SUCCESS;
    }
};

//! @}
}
