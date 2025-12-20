// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MathLib/LinAlg/MatrixSpecifications.h"

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
 */
class EquationSystem
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

    /*! Check whether normalization of A and rhs is required.
     *
     * \remark
     * In some processes, a normalization operation is required, to calculate
     * A^T * A, and overwrite A; also calculate A^T * rhs and overwrite rhs.
     * This parameter reflect whether such operation is required.
     */
    virtual bool requiresNormalization() const = 0;

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

    virtual MathLib::MatrixSpecifications getMatrixSpecifications(
        const int process_id) const = 0;

    virtual ~EquationSystem() = default;
};

//! @}
}  // namespace NumLib
