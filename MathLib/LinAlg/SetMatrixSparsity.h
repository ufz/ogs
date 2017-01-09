/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifndef MATHLIB_SETMATRIXSPARSITY_H_
#define MATHLIB_SETMATRIXSPARSITY_H_

namespace MathLib
{

/// Default implementation of SetMatrixSparsity class called by
/// setMatrixSparsity.
/// This is a workaround for partial function specialization.
template <typename MATRIX, typename SPARSITY_PATTERN>
struct SetMatrixSparsity
{
    void operator()(MATRIX&, SPARSITY_PATTERN const&)
    { }
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

} // MathLib

#ifdef USE_LIS
#include "Lis/LisMatrix.h"
#endif  // USE_LIS

#ifdef OGS_USE_EIGEN
#include "Eigen/EigenMatrix.h"
#endif  // OGS_USE_EIGEN

#endif  // MATHLIB_SETMATRIXSPARSITY_H_
