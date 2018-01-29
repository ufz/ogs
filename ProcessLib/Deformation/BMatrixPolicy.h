/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

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
