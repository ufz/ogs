// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cmath>
#include <optional>

#include "NumLib/Fem/AverageGradShapeFunction.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"

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
}  // namespace detail

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

    BMatrixType B = BMatrixType::Zero(
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim),
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

template <int DisplacementDim, int NPOINTS, typename ShapeFunction,
          typename BBarMatrixType, typename ShapeMatricesType, typename IpData>
BBarMatrixType computeDilatationalBbar(
    std::vector<IpData, Eigen::aligned_allocator<IpData>> const& ip_data,
    MeshLib::Element const& element,
    NumLib::GenericIntegrationMethod const& integration_method,
    const bool is_axially_symmetric)
{
    unsigned const n_integration_points =
        integration_method.getNumberOfPoints();

    BBarMatrixType B_bar = BBarMatrixType::Zero(3, ShapeFunction::NPOINTS);

    // Compute the element volume
    double volume = 0.0;
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& w = ip_data[ip].integration_weight;
        volume += w;
    }

    for (int i = 0; i < NPOINTS; i++)
    {
        B_bar.col(i).noalias() +=
            NumLib::averageGradShapeFunction<DisplacementDim, ShapeFunction,
                                             ShapeMatricesType, IpData>(
                i, element, integration_method, ip_data, is_axially_symmetric);
    }

    return B_bar / volume;
}

namespace detail
{
template <int DisplacementDim, int NPOINTS, typename BBarMatrixType,
          typename BMatrixType>
void applyBbar(BBarMatrixType const& B_bar, BMatrixType& B,
               const bool is_axially_symmetric)
{
    if constexpr (DisplacementDim == 3)
    {
        for (int i = 0; i < NPOINTS; ++i)
        {
            auto const B_bar_i = B_bar.col(i);

            // The following loop is based on the following facts:
            // B1 (dN/dx) is the 1st element of the 1st column of B,
            // B2 (dN/dy) is the 2nd element of the 2nd column of B,
            // B3 (dN/dz) is the 3rd element of the 3rd column of B.
            for (int k = 0; k < 3; k++)
            {
                // k-th column of B matrix at node i
                auto B_i_k = B.col(k * NPOINTS + i);

                B_i_k.template segment<3>(0) -=
                    Eigen::Vector3d::Constant((B_i_k[k] - B_bar_i[k]) / 3.0);
            }
        }

        return;
    }

    // 2D or axisymmetry
    for (int i = 0; i < NPOINTS; ++i)
    {
        auto B_i_0 = B.col(i);
        auto B_i_1 = B.col(NPOINTS + i);

        auto const B_bar_i = B_bar.col(i);

        if (is_axially_symmetric)
        {
            double const b0_dil_pertubation =
                (B_i_0[0] - B_bar_i[0] + B_i_0[2] - B_bar_i[2]);
            B_i_0.template segment<3>(0) -=
                Eigen::Vector3d::Constant((b0_dil_pertubation) / 3.);
            B_i_1.template segment<3>(0) -=
                Eigen::Vector3d::Constant((B_i_1[1] - B_bar_i[1]) / 3.);
            continue;
        }

        // Plane strain
        B_i_0.template segment<2>(0) -=
            Eigen::Vector2d::Constant((B_i_0[0] - B_bar_i[0]) / 2.);
        B_i_1.template segment<2>(0) -=
            Eigen::Vector2d::Constant((B_i_1[1] - B_bar_i[1]) / 2.);
    }
}
}  // namespace detail

/// Fills a B matrix, or a B bar matrix if required.
template <int DisplacementDim, int NPOINTS, typename BBarMatrixType,
          typename BMatrixType, typename N_Type, typename DNDX_Type>
BMatrixType computeBMatrixPossiblyWithBbar(
    DNDX_Type const& dNdx,
    N_Type const& N,
    std::optional<BBarMatrixType> const& B_dil_bar,
    const double radius,
    const bool is_axially_symmetric)
{
    auto B = computeBMatrix<DisplacementDim, NPOINTS, BMatrixType, N_Type,
                            DNDX_Type>(dNdx, N, radius, is_axially_symmetric);

    if (!B_dil_bar)
    {
        return B;
    }

    detail::applyBbar<DisplacementDim, NPOINTS, BBarMatrixType, BMatrixType>(
        *B_dil_bar, B, is_axially_symmetric);

    return B;
}

}  // namespace LinearBMatrix
}  // namespace ProcessLib
