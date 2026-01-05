// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace MathLib
{

/// Default implementation of SetMatrixSparsity class called by
/// setMatrixSparsity.
/// This is a workaround for partial function specialization.
template <typename MATRIX, typename SPARSITY_PATTERN>
struct SetMatrixSparsity
{
    void operator()(MATRIX& /*unused*/,
                    SPARSITY_PATTERN const& /*unused*/) const
    {
    }
};

/// Sets the sparsity pattern of the underlying matrix.
/// To allow partial specialization a SetMatrixSparsity template is
/// instantiated, to which the matrix and the sparsity_pattern are passed.
template <typename MATRIX, typename SPARSITY_PATTERN>
void setMatrixSparsity(MATRIX& matrix, SPARSITY_PATTERN const& sparsity_pattern)
{
    SetMatrixSparsity<MATRIX, SPARSITY_PATTERN> set_sparsity;
    set_sparsity(matrix, sparsity_pattern);
}

}  // namespace MathLib

#ifdef USE_LIS
#include "Lis/LisMatrix.h"
#endif  // USE_LIS

#include "Eigen/EigenMatrix.h"
