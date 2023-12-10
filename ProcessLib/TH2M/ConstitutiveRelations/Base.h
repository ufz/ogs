/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "MathLib/KelvinVector.h"

namespace ProcessLib::TH2M
{
template <int DisplacementDim>
using KelvinMatrix = MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

/// Convenience alias for not a number.
static constexpr double nan = std::numeric_limits<double>::quiet_NaN();

/// Used to set a Kelvin matrix to all not-a-number.
template <int DisplacementDim>
constexpr KelvinMatrix<DisplacementDim> KMnan()
{
    return KelvinMatrix<DisplacementDim>::Constant(nan);
}
}  // namespace ProcessLib::TH2M
