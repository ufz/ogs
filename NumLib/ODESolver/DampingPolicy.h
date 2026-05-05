// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

namespace NumLib
{
/// Optional plug-in for FixedDampingStrategy: given the current iterate and
/// the proposed Newton update, return a (possibly reduced) damping factor.
/// A null pointer means "no extra clamping"; the strategy uses its base
/// damping unchanged.
class DampingPolicy
{
public:
    virtual double apply(GlobalVector const& minus_delta_x,
                         GlobalVector const& x,
                         double const base_damping) const = 0;
    virtual ~DampingPolicy() = default;
};

}  // namespace NumLib
