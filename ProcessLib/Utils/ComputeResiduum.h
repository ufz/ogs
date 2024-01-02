/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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
