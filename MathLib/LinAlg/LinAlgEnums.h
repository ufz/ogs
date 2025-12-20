// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string>

namespace MathLib
{

/// Norm type. Not declared as class type in order to use the members as
/// integers.
enum class VecNormType
{
    NORM1,       ///< \f$\sum_i |x_i|\f$
    NORM2,       ///< \f$\sqrt(\sum_i (x_i)^2)\f$
    INFINITY_N,  ///< \f$\mathrm{max}_i |x_i|\f$
    INVALID
};

/// convert VecNormType to string
std::string convertVecNormTypeToString(VecNormType normType);

/// convert string to VecNormType
VecNormType convertStringToVecNormType(const std::string& str);

/// Determines if Dirichlet BCs will be applied properly to both the global
/// matrix A and the global rhs vector.
///
/// If the linear solver can reuse A from an earlier call, the Dirichlet BC
/// application to A can be incomplete, thereby faster.
enum class DirichletBCApplicationMode
{
    COMPLETE_MATRIX_UPDATE,        ///< Both A and b fully updated
    FAST_INCOMPLETE_MATRIX_UPDATE  ///< A partially updated, b fully updated
};

}  // end namespace MathLib
