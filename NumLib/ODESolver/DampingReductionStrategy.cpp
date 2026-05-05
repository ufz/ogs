// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "DampingReductionStrategy.h"

#include <algorithm>

#include "DampingPolicy.h"
#include "MathLib/LinAlg/LinAlg.h"

namespace NumLib
{
DampingReductionStrategy::DampingReductionStrategy(double damping,
                                                   double damping_reduction)
    : _damping(damping), _damping_reduction(damping_reduction)
{
}

StepResult DampingReductionStrategy::applyStep(
    GlobalVector const& x,
    GlobalVector const& minus_delta_x,
    GlobalVector const& /*res*/,
    GlobalMatrix const& /*J*/,
    GlobalVector& x_new,
    NewtonStepContext& /*ctx*/,
    int iteration)
{
    double damping = _damping;

    if (_damping_policy)
    {
        damping = _damping_policy->apply(minus_delta_x, x, _damping);
    }

    damping +=
        (1.0 - damping) * std::clamp(iteration / _damping_reduction, 0.0, 1.0);

    MathLib::LinAlg::axpy(x_new, -damping, minus_delta_x);

    return {.success = true, .step_length = damping, .x_new_is_set = true};
}

}  // namespace NumLib
