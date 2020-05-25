/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MathLib/KelvinVector.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

namespace ProcessLib
{
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
    static int const number_of_dof_ = ShapeFunction::NPOINTS * DisplacementDim;
    static int const kelvin_vector_size_ =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;

public:
    using StiffnessMatrixType = MatrixType<number_of_dof_, number_of_dof_>;

    /// Rhs residual
    using NodalForceVectorType = VectorType<number_of_dof_>;

    /// This type can be different (fixed vs. dynamic size) from the
    /// MathLib::KelvinVector::KelvinVectorType.
    using KelvinVectorType = VectorType<kelvin_vector_size_>;

    /// This type can be different (fixed vs. dynamic size) from the
    /// MathLib::KelvinVector::KelvinMatrixType.
    using KelvinMatrixType =
        MatrixType<kelvin_vector_size_, kelvin_vector_size_>;

    using BMatrixType = MatrixType<kelvin_vector_size_, number_of_dof_>;
};
}  // namespace ProcessLib
