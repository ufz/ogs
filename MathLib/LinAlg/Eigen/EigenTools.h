/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "BaseLib/Error.h"
#include "EigenMatrix.h"  // for EigenMatrix::IndexType

namespace MathLib
{
class EigenVector;

/**
 * apply known solutions to a system of linear equations
 *
 * @param A                 Coefficient matrix
 * @param b                 RHS vector
 * @param vec_knownX_id    a vector of known solution entry IDs
 * @param vec_knownX_x     a vector of known solutions
 * @param penalty_scaling value for scaling some matrix and right hand side
 * entries to enforce some conditions, value ignored in the current
 * implementation
 */
void applyKnownSolution(
    EigenMatrix& A, EigenVector& b, EigenVector& /*x*/,
    const std::vector<EigenMatrix::IndexType>& vec_knownX_id,
    const std::vector<double>& vec_knownX_x, double penalty_scaling = 1e+10);

}  // namespace MathLib
