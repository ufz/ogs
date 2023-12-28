/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "MathLib/VectorizedTensor.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

namespace ProcessLib
{
/// An implementation of G-Matrix policy using same matrix and vector types
/// (fixed size or dynamic) as in the ShapeMatrixPolicyType.
template <typename ShapeFunction, int DisplacementDim>
class GMatrixPolicyType
{
private:
    /// Fixed size vector type independent of the ShapeMatrixPolicy needed for
    /// storage of gradient vector in the MPL VariableArray.
    template <int N>
    using VectorTypeFixedSize = typename EigenFixedShapeMatrixPolicy<
        ShapeFunction, DisplacementDim>::template VectorType<N>;

    /// Reusing the ShapeMatrixPolicy matrix type.
    template <int N, int M>
    using MatrixType = typename ShapeMatrixPolicyType<
        ShapeFunction, DisplacementDim>::template MatrixType<N, M>;

    // Dimensions of specific b-matrix for n-points and displacement dimension.
    static int const _number_of_dof = ShapeFunction::NPOINTS * DisplacementDim;

public:
    // For the 2D case the 33-component is needed (and the four entries
    // of the non-symmetric matrix); In 3d there are nine entries.
    using GradientMatrixType =
        MatrixType<MathLib::VectorizedTensor::size(DisplacementDim),
                   _number_of_dof>;
    using GradientVectorType =
        VectorTypeFixedSize<MathLib::VectorizedTensor::size(DisplacementDim)>;
};
}  // namespace ProcessLib
