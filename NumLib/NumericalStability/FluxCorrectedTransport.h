/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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
