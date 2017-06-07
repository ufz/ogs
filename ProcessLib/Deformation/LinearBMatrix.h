/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cmath>

namespace ProcessLib
{
namespace LinearBMatrix
{
namespace detail
{
template <int NPOINTS, typename DNDX_Type, typename BMatrixType>
void fillBMatrix2DCartesianPart(DNDX_Type const& dNdx, BMatrixType& B)
{
    for (int i = 0; i < NPOINTS; ++i)
    {
        B(1, NPOINTS + i) = dNdx(1, i);
        B(3, i) = dNdx(1, i) / std::sqrt(2);
        B(3, NPOINTS + i) = dNdx(0, i) / std::sqrt(2);
        B(0, i) = dNdx(0, i);
    }
}
}  // detail

/// Fills a B-matrix based on given shape function dN/dx values.
template <int DisplacementDim,
          int NPOINTS,
          typename BMatrixType,
          typename N_Type,
          typename DNDX_Type>
BMatrixType computeBMatrix(DNDX_Type const& dNdx,
                           N_Type const& N,
                           const double radius,
                           const bool is_axially_symmetric)
{
    static_assert(0 < DisplacementDim && DisplacementDim <= 3,
                  "LinearBMatrix::computeBMatrix: DisplacementDim must be in "
                  "range [1,3].");

    BMatrixType B =
        BMatrixType::Zero(KelvinVectorDimensions<DisplacementDim>::value,
                          NPOINTS * DisplacementDim);

    switch (DisplacementDim)
    {
        case 3:
            for (int i = 0; i < NPOINTS; ++i)
            {
                B(2, 2 * NPOINTS + i) = dNdx(2, i);
                B(4, NPOINTS + i) = dNdx(2, i) / std::sqrt(2);
                B(4, 2 * NPOINTS + i) = dNdx(1, i) / std::sqrt(2);
                B(5, i) = dNdx(2, i) / std::sqrt(2);
                B(5, 2 * NPOINTS + i) = dNdx(0, i) / std::sqrt(2);
            }
            detail::fillBMatrix2DCartesianPart<NPOINTS>(dNdx, B);
            break;
        case 2:
            detail::fillBMatrix2DCartesianPart<NPOINTS>(dNdx, B);
            if (is_axially_symmetric)
            {
                for (int i = 0; i < NPOINTS; ++i)
                {
                    B(2, i) = N[i] / radius;
                }
            }
            break;
        default:
            break;
    }

    return B;
}

}  // namespace LinearBMatrix
}  // namespace ProcessLib
