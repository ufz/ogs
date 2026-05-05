// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "NewtonStepStrategy.h"

namespace NumLib
{
/*! Newton step strategy that applies a fixed base damping factor and relaxes
 *  it linearly toward 1.0 as iterations proceed.
 *
 * The effective damping at iteration \f$ k \f$ is
 * \f[
 *   \alpha_k = \alpha + (1 - \alpha)\,\mathrm{clamp}(k / r,\, 0,\, 1),
 * \f]
 * where \f$ \alpha \f$ is the base damping factor and \f$ r \f$ is the
 * reduction parameter. The accepted step is
 * \f$ x_{\rm new} = x - \alpha_k \cdot \Delta x \f$.
 */
class DampingReductionStrategy final : public NewtonStepStrategy
{
public:
    DampingReductionStrategy(double damping, double damping_reduction);

    StepResult applyStep(GlobalVector const& x,
                         GlobalVector const& minus_delta_x,
                         GlobalVector const& res,
                         GlobalMatrix const& J,
                         GlobalVector& x_new,
                         NewtonStepContext& ctx,
                         int iteration) override;

private:
    double _damping;
    double _damping_reduction;
};

}  // namespace NumLib
