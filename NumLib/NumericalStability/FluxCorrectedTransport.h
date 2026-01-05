// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MathLib/LinAlg/MatrixSpecifications.h"
#include "MathLib/LinAlg/MatrixVectorTraits.h"
#include "NumericalStabilization.h"

namespace NumLib
{
void computeFluxCorrectedTransport(
    NumericalStabilization const& stabilizer, const double t, const double dt,
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    const MathLib::MatrixSpecifications& matrix_specification, GlobalMatrix& M,
    GlobalMatrix& K, GlobalVector& b);
}  // namespace NumLib
