/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cmath>

#include "ProcessLib/Deformation/BMatrixPolicy.h"

namespace ProcessLib
{
namespace NonLinearBMatrix
{
/// Fills a non linear B-matrix based on given shape function dN/dx values and
/// displacement gradient.
template <int DisplacementDim,
          int NPOINTS,
          typename BMatrixType,
          typename GradientVectorType,
          typename N_Type,
          typename DNDX_Type>
BMatrixType computeBMatrix(DNDX_Type const& dNdx,
                           N_Type const& N,
                           GradientVectorType const& grad_u,
                           const double radius,
                           const bool is_axially_symmetric)
{
    static_assert(0 < DisplacementDim && DisplacementDim <= 3,
                  "NonLinearBMatrix::computeBMatrix: DisplacementDim must be "
                  "in range [1,3].");

    BMatrixType B = BMatrixType::Zero(
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim),
        NPOINTS * DisplacementDim);

    switch (DisplacementDim)
    {
        case 3:
            for (int i = 0; i < NPOINTS; ++i)
            {
                B(0, i) = dNdx(0, i) * grad_u[0];
                B(1, i) = dNdx(1, i) * grad_u[1];
                B(2, i) = dNdx(2, i) * grad_u[2];

                B(0, NPOINTS + i) = dNdx(0, i) * grad_u[3];
                B(1, NPOINTS + i) = dNdx(1, i) * grad_u[4];
                B(2, NPOINTS + i) = dNdx(2, i) * grad_u[5];

                B(0, 2 * NPOINTS + i) = dNdx(0, i) * grad_u[6];
                B(1, 2 * NPOINTS + i) = dNdx(1, i) * grad_u[7];
                B(2, 2 * NPOINTS + i) = dNdx(2, i) * grad_u[8];

                B(3, i) = (dNdx(1, i) * grad_u[0] + dNdx(0, i) * grad_u[1]) /
                          std::sqrt(2);
                B(4, i) = (dNdx(2, i) * grad_u[1] + dNdx(1, i) * grad_u[2]) /
                          std::sqrt(2);
                B(5, i) = (dNdx(2, i) * grad_u[0] + dNdx(0, i) * grad_u[2]) /
                          std::sqrt(2);

                B(3, NPOINTS + i) =
                    (dNdx(1, i) * grad_u[3] + dNdx(0, i) * grad_u[4]) /
                    std::sqrt(2);
                B(4, NPOINTS + i) =
                    (dNdx(2, i) * grad_u[4] + dNdx(1, i) * grad_u[5]) /
                    std::sqrt(2);
                B(5, NPOINTS + i) =
                    (dNdx(2, i) * grad_u[3] + dNdx(0, i) * grad_u[5]) /
                    std::sqrt(2);

                B(3, 2 * NPOINTS + i) =
                    (dNdx(1, i) * grad_u[6] + dNdx(0, i) * grad_u[7]) /
                    std::sqrt(2);
                B(4, 2 * NPOINTS + i) =
                    (dNdx(2, i) * grad_u[7] + dNdx(1, i) * grad_u[8]) /
                    std::sqrt(2);
                B(5, 2 * NPOINTS + i) =
                    (dNdx(2, i) * grad_u[6] + dNdx(0, i) * grad_u[8]) /
                    std::sqrt(2);
            }
            break;
        case 2:
            for (int i = 0; i < NPOINTS; ++i)
            {
                B(0, i) = dNdx(0, i) * grad_u[0];
                B(1, i) = dNdx(1, i) * grad_u[1];
                // B(2, i) = 0;

                B(0, NPOINTS + i) = dNdx(0, i) * grad_u[2];
                B(1, NPOINTS + i) = dNdx(1, i) * grad_u[3];
                // B(2, NPOINTS + i) = 0;

                B(3, i) = (dNdx(1, i) * grad_u[0] + dNdx(0, i) * grad_u[1]) /
                          std::sqrt(2);
                B(3, NPOINTS + i) =
                    (dNdx(1, i) * grad_u[2] + dNdx(0, i) * grad_u[3]) /
                    std::sqrt(2);
            }
            if (is_axially_symmetric)
            {
                for (int i = 0; i < NPOINTS; ++i)
                {
                    B(2, i) = grad_u[4] * N[i] / radius;
                }
            }
            break;
        default:
            break;
    }

    return B;
}
}  // namespace NonLinearBMatrix
}  // namespace ProcessLib
