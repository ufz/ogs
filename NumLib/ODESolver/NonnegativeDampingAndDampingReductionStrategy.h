// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <optional>

#include "NewtonStepStrategy.h"

namespace NumLib
{
/*! Newton step strategy that applies a fixed damping factor to every update.
 *
 * The accepted step is \f$ x_{\rm new} = x - \alpha \cdot \Delta x \f$, where
 * \f$ \alpha \f$ is the damping factor (1.0 = full Newton step, no damping).
 * An optional \c damping_reduction parameter can be used to relax the damping
 * linearly towards 1.0 as iterations proceed.
 */
class FixedDampingStrategy final : public NewtonStepStrategy
{
public:
    /*! Constructs a fixed-damping step strategy.
     *
     * \param damping          Initial damping factor; must be positive.
     *                         1.0 gives an undamped Newton step.
     * \param damping_reduction If set, the effective damping is increased
     *                         towards 1.0 by \c iteration / damping_reduction
     *                         each iteration (clamped to [damping, 1]).
     */
    FixedDampingStrategy(
        double damping, std::optional<double> damping_reduction);

    StepResult applyStep(GlobalVector const& x,
                         GlobalVector const& minus_delta_x,
                         GlobalVector const& res,
                         GlobalMatrix const& J,
                         GlobalVector& x_new,
                         NewtonStepContext& ctx,
                         int iteration) override;

    /// Injects the convergence criterion so that applyStep() can query
    /// hasNonNegativeDamping() and getDampingFactor().
    void setConvergenceCriterion(ConvergenceCriterion& conv_crit) override
    {
        _conv_crit = &conv_crit;
    }

private:
    double _damping;  //!< Base damping factor.
    std::optional<double>
        _damping_reduction;  //!< Controls per-iteration damping relaxation.
    /// Convergence criterion used to query getDampingFactor() and
    /// hasNonNegativeDamping(). Set via setConvergenceCriterion().
    ConvergenceCriterion* _conv_crit = nullptr;
};

}  // namespace NumLib
