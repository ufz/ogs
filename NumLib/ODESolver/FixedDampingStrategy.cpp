// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "FixedDampingStrategy.h"

#include "DampingPolicy.h"
#include "MathLib/LinAlg/LinAlg.h"

namespace NumLib
{
FixedDampingStrategy::FixedDampingStrategy(double damping) : _damping(damping)
{
}

StepResult FixedDampingStrategy::applyStep(GlobalVector const& x,
                                           GlobalVector const& minus_delta_x,
                                           GlobalVector const& /*res*/,
                                           GlobalMatrix const& /*J*/,
                                           GlobalVector& x_new,
                                           NewtonStepContext& /*ctx*/,
                                           int /*iteration*/)
{
    double damping = _damping;

    if (_damping_policy)
    {
        damping = _damping_policy->apply(minus_delta_x, x, _damping);
    }

    MathLib::LinAlg::axpy(x_new, -damping, minus_delta_x);

    return {.success = true, .step_length = damping, .x_new_is_set = true};
}

}  // namespace NumLib
