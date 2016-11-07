/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_LIE_COMMON_HMATRIXPOLICYTYPE_H_
#define PROCESSLIB_LIE_COMMON_HMATRIXPOLICYTYPE_H_

#include "NumLib/Fem/ShapeMatrixPolicy.h"

namespace ProcessLib
{

/// An implementation of H-Matrix policy using same matrix and vector types
/// (fixed size or dynamic) as in the ShapeMatrixPolicyType.
template <typename ShapeFunction, unsigned DisplacementDim>
class HMatrixPolicyType
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

    // Dimensions of specific H-matrix for n-points and displacement dimension.
    static int const _number_of_dof = ShapeFunction::NPOINTS * DisplacementDim;

public:
    using StiffnessMatrixType = MatrixType<_number_of_dof, _number_of_dof>;
    using NodalForceVectorType = VectorType<_number_of_dof>;

    using HMatrixType = MatrixType<DisplacementDim, _number_of_dof>;

    using ConstitutiveMatrixType = MatrixType<DisplacementDim, DisplacementDim>;
    using ForceVectorType = VectorType<DisplacementDim>;
};


/// Fills a H-matrix based on given shape function
template <int DisplacementDim, int NPOINTS, typename N_Type,
          typename HMatrixType>
void computeHMatrix(N_Type const& N, HMatrixType& H)
{
    static_assert(1 < DisplacementDim && DisplacementDim <= 3,
                  "LinearHMatrix::computeHMatrix: DisplacementDim must be in "
                  "range (1,3].");

    H.setZero();

    for (unsigned j=0; j<DisplacementDim; j++)
        H.block(j, j*NPOINTS, 1, NPOINTS) = N;
}

}  // namespace ProcessLib

#endif // PROCESSLIB_LIE_COMMON_HMATRIXPOLICYTYPE_H_
