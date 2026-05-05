// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "NewtonStepStrategy.h"

namespace NumLib
{
/*! Newton step strategy that applies a fixed damping factor to every update.
 *
 * The accepted step is \f$ x_{\rm new} = x - \alpha \cdot \Delta x \f$, where
 * \f$ \alpha \f$ is the damping factor (1.0 = full Newton step, no damping).
 */
class FixedDampingStrategy final : public NewtonStepStrategy
{
public:
    explicit FixedDampingStrategy(double damping);

    StepResult applyStep(GlobalVector const& x,
                         GlobalVector const& minus_delta_x,
                         GlobalVector const& res,
                         GlobalMatrix const& J,
                         GlobalVector& x_new,
                         NewtonStepContext& ctx,
                         int iteration) override;

private:
    double _damping;
};

}  // namespace NumLib
