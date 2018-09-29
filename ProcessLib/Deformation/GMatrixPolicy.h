/**
 * \file
 *
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
/// An implementation of G-Matrix policy using same matrix and vector types
/// (fixed size or dynamic) as in the ShapeMatrixPolicyType.
template <typename ShapeFunction, unsigned DisplacementDim>
class GMatrixPolicyType
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

public:
    // For the 2D case the 33-component is needed (and the four entries
    // of the non-symmetric matrix); In 3d there are nine entries.
    using GradientMatrixType = MatrixType<DisplacementDim * DisplacementDim +
                                              (DisplacementDim == 2 ? 1 : 0),
                                          _number_of_dof>;
    using GradientVectorType = VectorType<DisplacementDim * DisplacementDim +
                                          (DisplacementDim == 2 ? 1 : 0)>;
};
}  // namespace ProcessLib
