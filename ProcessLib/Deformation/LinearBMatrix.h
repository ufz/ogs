/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_FEM_DEFORMATION_LINEARBMATRIX_H_
#define PROCESSLIB_FEM_DEFORMATION_LINEARBMATRIX_H_

#include <cmath>

namespace ProcessLib
{
namespace LinearBMatrix
{
namespace detail
{
template <int NPOINTS, typename DNDX_Type, typename BMatrixType>
void fillBMatrix2DCartesianPart(DNDX_Type const& dNdx, BMatrixType& b_matrix)
{
    for (int i = 0; i < NPOINTS; ++i) {
        b_matrix(1, NPOINTS + i) = dNdx(1, i);
        b_matrix(3, i) = dNdx(1, i) / std::sqrt(2);
        b_matrix(3, NPOINTS + i) = dNdx(0, i) / std::sqrt(2);
        b_matrix(0, i) = dNdx(0, i);
    }
}
}  // detail

/// Fills a B-matrix based on given shape function dN/dx values.
template <int DisplacementDim,
          int NPOINTS,
          typename N_Type,
          typename DNDX_Type,
          typename BMatrixType>
void computeBMatrix(DNDX_Type const& dNdx,
                    BMatrixType& b_matrix,
                    const bool is_axially_symmetric,
                    N_Type const& N,
                    const double radius)
{
    static_assert(0 < DisplacementDim && DisplacementDim <= 3,
                  "LinearBMatrix::computeBMatrix: DisplacementDim must be in "
                  "range [1,3].");

    b_matrix.setZero();

    switch (DisplacementDim)
    {
        case 3:
            for (int i = 0; i < NPOINTS; ++i)
            {
                b_matrix(2, 2 * NPOINTS + i) = dNdx(2, i);
                b_matrix(4, NPOINTS + i) = dNdx(2, i) / std::sqrt(2);
                b_matrix(4, 2 * NPOINTS + i) = dNdx(1, i) / std::sqrt(2);
                b_matrix(5, i) = dNdx(2, i) / std::sqrt(2);
                b_matrix(5, 2 * NPOINTS + i) = dNdx(0, i) / std::sqrt(2);
            }
            detail::fillBMatrix2DCartesianPart<NPOINTS>(dNdx, b_matrix);
            break;
        case 2:
            detail::fillBMatrix2DCartesianPart<NPOINTS>(dNdx, b_matrix);
            if (is_axially_symmetric)
            {
                for (int i = 0; i < NPOINTS; ++i)
                {
                    b_matrix(2, i) = N[i] / radius;
                }
            }
            break;
        default:
            break;
    }
}

}  // namespace LinearBMatrix
}  // namespace ProcessLib

#endif  // PROCESSLIB_FEM_DEFORMATION_LINEARBMATRIX_H_
