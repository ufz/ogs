// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "FixedDampingStrategy.h"

#include "DampingPolicy.h"

namespace NumLib
{
FixedDampingStrategy::FixedDampingStrategy(double damping) : _damping(damping)
{
}

StepResult FixedDampingStrategy::applyStep(GlobalVector const& minus_delta_x,
                                           GlobalMatrix const& /*J*/,
                                           GlobalVector& x_new,
                                           NewtonStepContext const& ctx)
{
    double const damping =
        applyDampingPolicy(minus_delta_x, *ctx.x[ctx.process_id], _damping);

    applyDampedStep(x_new, minus_delta_x, damping);

    return {.success = true, .step_length = damping, .x_new_is_set = true};
}

}  // namespace NumLib
