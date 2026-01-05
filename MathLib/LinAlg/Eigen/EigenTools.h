// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

#include "BaseLib/Error.h"
#include "EigenMatrix.h"  // for EigenMatrix::IndexType
#include "MathLib/LinAlg/LinAlgEnums.h"

namespace MathLib
{
class EigenVector;

/**
 * apply known solutions to a system of linear equations
 *
 * \param A                 Coefficient matrix
 * \param b                 RHS vector
 * \param vec_knownX_id    a vector of known solution entry IDs
 * \param vec_knownX_x     a vector of known solutions
 * \param mode             determines if the modification of matrix A should be
 * "complete", incomplete is faster but leaves A in a garbage state and only
 * modifies b properly
 */
void applyKnownSolution(
    EigenMatrix& A, EigenVector& b, EigenVector& /*x*/,
    const std::vector<EigenMatrix::IndexType>& vec_knownX_id,
    const std::vector<double>& vec_knownX_x,
    DirichletBCApplicationMode const mode);

}  // namespace MathLib
