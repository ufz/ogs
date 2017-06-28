/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include "MathLib/LinAlg/LinAlg.h" // For MathLib::VecNormType
#include "NumLib/NumericsConfig.h"

namespace BaseLib {
class ConfigTree;
}  // BaseLib

namespace NumLib
{
/*! Convergence criterion for iterative algorithms, like nonlinear solvers.
 *
 * It is able to check the difference of the solution vector between iterations
 * and the residual of the equation being solved (e.g. the nonlinear equation
 * system).
 */
class ConvergenceCriterion
{
public:
    ConvergenceCriterion(const MathLib::VecNormType norm_type)
        : _norm_type(norm_type)
    {
    }

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
    //! \param minus_delta_x the current solution update
    //! \param x the new solution from the current iteration
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
    virtual void preFirstIteration() { _is_first_iteration = true; }

    //! Tell the ConvergenceCriterion that it is not called for the first time
    //! (while solving a coupling system).
    virtual void setNoFirstIteration() { _is_first_iteration = false; }

    //! Indicate that a new iteration now starts.
    //!
    //! A concrete implementation of ConvergenceCriterion might want to check
    //! both delta x and the residual. I.e., in a new iteration both checks have
    //! to be done anew. This method will make the ConvergenceCriterion forget
    //! the result of all checks already done, s.t. all necessary checks will
    //! have to be repeated in order to satisfy the ConvergenceCriterion.
    virtual void reset() { _satisfied = true; _is_first_iteration = false; };

    //! Tell if the convergence criterion is satisfied.
    virtual bool isSatisfied() const { return _satisfied; };

    MathLib::VecNormType getVectorNormType () const { return _norm_type; }

    virtual ~ConvergenceCriterion() = default;

protected:
    bool _satisfied = true;
    bool _is_first_iteration = true;
    const MathLib::VecNormType _norm_type;
};

//! Creates a convergence criterion from the given configuration.
std::unique_ptr<ConvergenceCriterion> createConvergenceCriterion(
    BaseLib::ConfigTree const& config);

//! Returns if |numerator/denominator| < |reltol|.
//! This method copes with the case that denominator = 0 by always adding
//! epsilon to the denominator.
bool checkRelativeTolerance(double const reltol,
                            double const numerator,
                            double const denominator);

} // namespace NumLib
