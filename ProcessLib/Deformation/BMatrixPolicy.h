// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MathLib/KelvinVector.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

namespace ProcessLib
{
/// An implementation of B-Matrix policy using same matrix and vector types
/// (fixed size or dynamic) as in the ShapeMatrixPolicyType.
template <typename ShapeFunction, int DisplacementDim>
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
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);

public:
    using StiffnessMatrixType = MatrixType<_number_of_dof, _number_of_dof>;

    /// Rhs residual
    using NodalForceVectorType = VectorType<_number_of_dof>;

    /// This type can be different (fixed vs. dynamic size) from the
    /// MathLib::KelvinVector::KelvinVectorType.
    using KelvinVectorType = VectorType<_kelvin_vector_size>;

    /// This type can be different (fixed vs. dynamic size) from the
    /// MathLib::KelvinVector::KelvinMatrixType.
    using KelvinMatrixType =
        MatrixType<_kelvin_vector_size, _kelvin_vector_size>;

    using BMatrixType = MatrixType<_kelvin_vector_size, _number_of_dof>;

    using BBarMatrixType = MatrixType<3, ShapeFunction::NPOINTS>;
};
}  // namespace ProcessLib
