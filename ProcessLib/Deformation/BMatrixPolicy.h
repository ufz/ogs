/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/Fem/ShapeMatrixPolicy.h"

namespace ProcessLib
{
/// Kelvin vector dimensions for given displacement dimension.
template <int DisplacementDim>
struct KelvinVectorDimensions;

template <>
struct KelvinVectorDimensions<2>
{
    static int const value = 4;
};

template <>
struct KelvinVectorDimensions<3>
{
    static int const value = 6;
};

//
// Kelvin vector and matrix templates for given displacement dimension.
//

/// \todo Maybe better to move the KelvinVector/MatrixType into
/// MaterialLib/SolidModels/KelvinVector.h.

/// Kelvin vector type for given displacement dimension.
/// \note The Eigen vector is always a fixed size vector in contrast to the
/// BMatrixPolicyType::KelvinVectorType.
template <int DisplacementDim>
using KelvinVectorType =
    Eigen::Matrix<double, KelvinVectorDimensions<DisplacementDim>::value, 1,
                  Eigen::ColMajor>;

/// Kelvin matrix type for given displacement dimension.
/// \note The Eigen matrix is always a fixed size matrix in contrast to the
/// BMatrixPolicyType::KelvinMatrixType.
template <int DisplacementDim>
using KelvinMatrixType =
    Eigen::Matrix<double, KelvinVectorDimensions<DisplacementDim>::value,
                  KelvinVectorDimensions<DisplacementDim>::value,
                  Eigen::RowMajor>;

/// An implementation of B-Matrix policy using same matrix and vector types
/// (fixed size or dynamic) as in the ShapeMatrixPolicyType.
template <typename ShapeFunction, unsigned DisplacementDim>
class BMatrixPolicyType
{
private:
    /// Reusing the ShapeMatrixPolicy vector type.
    template <int N>
    using VectorType =
        typename ShapeMatrixPolicyType<ShapeFunction,
                                       DisplacementDim>::template VectorType<N>;

    /// Reusing the ShapeMatrixPolicy matrix type.
    template <int N, int M>
    using MatrixType = typename ShapeMatrixPolicyType<
        ShapeFunction, DisplacementDim>::template MatrixType<N, M>;

    // Dimensions of specific b-matrix for n-points and displacement dimension.
    static int const _number_of_dof = ShapeFunction::NPOINTS * DisplacementDim;
    static int const _kelvin_vector_size =
        KelvinVectorDimensions<DisplacementDim>::value;

public:
    using StiffnessMatrixType = MatrixType<_number_of_dof, _number_of_dof>;

    /// Rhs residual
    using NodalForceVectorType = VectorType<_number_of_dof>;

    /// This type can be different (fixed vs. dynamic size) from the
    /// ProcessLib::KelvinVectorType.
    using KelvinVectorType = VectorType<_kelvin_vector_size>;

    /// This type can be different (fixed vs. dynamic size) from the
    /// ProcessLib::KelvinMatrixType.
    using KelvinMatrixType =
        MatrixType<_kelvin_vector_size, _kelvin_vector_size>;

    using BMatrixType = MatrixType<_kelvin_vector_size, _number_of_dof>;
};
}  // namespace ProcessLib
