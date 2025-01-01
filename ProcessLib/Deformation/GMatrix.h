/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <cmath>

namespace ProcessLib
{
namespace Deformation
{
/// Fills a G-matrix based on given shape function dN/dx values.
template <int DisplacementDim,
          int NPOINTS,
          typename N_Type,
          typename DNDX_Type,
          typename GMatrixType>
void computeGMatrix(DNDX_Type const& dNdx,
                    GMatrixType& g_matrix,
                    const bool is_axially_symmetric,
                    N_Type const& N,
                    const double radius)
{
    static_assert(0 < DisplacementDim && DisplacementDim <= 3,
                  "LinearGMatrix::computeGMatrix: DisplacementDim must be in "
                  "range [1,3].");

    g_matrix.setZero();

    switch (DisplacementDim)
    {
        case 3:
            // The gradient coordinates are organized in the following order:
            // (1,1), (1,2), (1,3)
            // (2,1), (2,2), (2,3)
            // (3,1), (3,2), (3,3)
            for (int i = 0; i < NPOINTS; ++i)
            {
                auto G_col0 = g_matrix.col(i);
                auto G_col1 = g_matrix.col(i + NPOINTS);
                auto G_col2 = g_matrix.col(i + 2 * NPOINTS);
                auto const dNidx = dNdx.col(i);

                G_col0.template segment<DisplacementDim>(0).noalias() = dNidx;
                G_col1.template segment<DisplacementDim>(DisplacementDim)
                    .noalias() = dNidx;
                G_col2.template segment<DisplacementDim>(2 * DisplacementDim)
                    .noalias() = dNidx;
            }

            break;
        case 2:
            // The gradient coordinates are organized in the following order:
            // (1,1), (1,2)
            // (2,1), (2,2)
            // (3,3)
            for (int i = 0; i < NPOINTS; ++i)
            {
                auto G_col0 = g_matrix.col(i);
                auto G_col1 = g_matrix.col(i + NPOINTS);
                auto const dNidx = dNdx.col(i);

                G_col0.template segment<DisplacementDim>(0).noalias() = dNidx;
                G_col1.template segment<DisplacementDim>(DisplacementDim)
                    .noalias() = dNidx;

                if (is_axially_symmetric)
                {
                    G_col0[4] = N[i] / radius;
                }
            }
            break;
        default:
            break;
    }
}

}  // namespace Deformation
}  // namespace ProcessLib
