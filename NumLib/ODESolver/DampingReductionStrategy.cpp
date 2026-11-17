// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "DampingReductionStrategy.h"

#include <algorithm>

#include "DampingPolicy.h"

namespace NumLib
{
DampingReductionStrategy::DampingReductionStrategy(double damping,
                                                   double damping_reduction)
    : _damping(damping), _damping_reduction(damping_reduction)
{
}

StepResult DampingReductionStrategy::applyStep(
    GlobalVector const& minus_delta_x,
    GlobalMatrix const& /*J*/,
    GlobalVector& x_new,
    NewtonStepContext const& ctx)
{
    double damping =
        applyDampingPolicy(minus_delta_x, *ctx.x[ctx.process_id], _damping);

    damping += (1.0 - damping) *
               std::clamp(ctx.iteration / _damping_reduction, 0.0, 1.0);

    applyDampedStep(x_new, minus_delta_x, damping);

    return {.success = true, .step_length = damping, .x_new_is_set = true};
}

}  // namespace NumLib
