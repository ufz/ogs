/**
 * \file
 *
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
            for (int d = 0; d < DisplacementDim; ++d)
            {
                for (int i = 0; i < NPOINTS; ++i)
                {
                    g_matrix(d + 0 * DisplacementDim, i + 0 * NPOINTS) =
                        dNdx(d, i);
                    g_matrix(d + 1 * DisplacementDim, i + 1 * NPOINTS) =
                        dNdx(d, i);
                    g_matrix(d + 2 * DisplacementDim, i + 2 * NPOINTS) =
                        dNdx(d, i);
                }
            }
            break;
        case 2:
            // The gradient coordinates are organized in the following order:
            // (1,1), (1,2)
            // (2,1), (2,2)
            // (3,3)
            for (int d = 0; d < DisplacementDim; ++d)
            {
                for (int i = 0; i < NPOINTS; ++i)
                {
                    g_matrix(d, i) = dNdx(d, i);
                    g_matrix(d + DisplacementDim, i + NPOINTS) = dNdx(d, i);
                }
            }
            if (is_axially_symmetric)
            {
                for (int i = 0; i < NPOINTS; ++i)
                {
                    g_matrix(4, i) = N[i] / radius;
                }
            }
            break;
        default:
            break;
    }
}

}  // namespace Deformation
}  // namespace ProcessLib
