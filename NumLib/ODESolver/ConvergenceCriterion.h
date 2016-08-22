/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_CONVERGENCECRITERION_H
#define NUMLIB_CONVERGENCECRITERION_H

#include <memory>
#include "NumLib/NumericsConfig.h"

namespace BaseLib {
class ConfigTree;
}  // BaseLib

namespace NumLib
{
/*! Convergence criterion for nonlinear solvers.
 *
 * It is able to check the difference of the solution vector between iterations
 * and the residual of the nonlinear equation system.
 */
class ConvergenceCriterion
{
public:
    //! Tells if the change of the solution between iterations is checked.
    //!
    //! \remark
    //! This method allows to save some computations if no such check will be
    //! done.
    virtual bool hasDeltaXCheck() const = 0;

    //! Tells if the residual is checked.
    //!
    //! \remark
    //! This method allows to save some computations if no such check will be
    //! done.
    virtual bool hasResidualCheck() const = 0;

    //! Check if the change of the solution between iterations satisfies the
    //! convergence criterion.
    //!
    //! \remark
    //! The Newton-Raphson solver computes \c minus_delta_x. \c x is needed for
    //! relative tolerances.
    virtual void checkDeltaX(GlobalVector const& minus_delta_x,
                             GlobalVector const& x) = 0;

    //! Check if the residual satisfies the convergence criterion.
    virtual void checkResidual(GlobalVector const& residual) = 0;

    //! Tell the ConvergenceCriterion that it is called for the first time now
    //! (while solving a specific nonlinear system).
    virtual void preFirstIteration() {}

    //! Indicate that a new iteration now starts.
    //!
    //! A concrete implementation of ConvergenceCriterion might want to check
    //! both delta x and the residual. I.e., in a new iteration both checks have
    //! to be done anew. This method will make the ConvergenceCriterion forget
    //! the result of all checks already done, s.t. all necessary checks will
    //! have to be repeated in order to satisfy the ConvergenceCriterion.
    virtual void reset() = 0;

    //! Tell if the convergence criterion is satisfied.
    virtual bool isSatisfied() const = 0;

    ~ConvergenceCriterion() = default;
};

//! Creates a convergence criterion from the given configuration.
std::unique_ptr<ConvergenceCriterion> createConvergenceCriterion(
    BaseLib::ConfigTree const& config);

} // namespace NumLib

#endif  // NUMLIB_CONVERGENCECRITERION_H
