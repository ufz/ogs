// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

namespace ProcessLib
{
/// Computes the residuum r = -M*x_dot - K*x + b. Negation for consistency with
/// the Newton scheme.
GlobalVector computeResiduum(double const dt, GlobalVector const& x,
                             GlobalVector const& x_prev, GlobalMatrix const& M,
                             GlobalMatrix const& K, GlobalVector const& b);

}  // namespace ProcessLib
