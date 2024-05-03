/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on May 3, 2024, 10:47 AM
 */

#pragma once

#include <cmath>
#include <tuple>

#include "GMatrix.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/VectorizedTensor.h"
#include "MeshLib/Elements/Element.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/Integration/GenericIntegrationMethod.h"

namespace ProcessLib
{
namespace NonLinearFbar
{
template <int DisplacementDim, typename GradientVectorType,
          typename GradientMatrixType, typename GlobalDimNodalMatrixType,
          typename ShapeFunction, typename ShapeMatricesType>
std::tuple<double, GlobalDimNodalMatrixType> computeFBarInitialVariables(
    bool const compute_detF0_only, Eigen::VectorXd const& u,
    NumLib::GenericIntegrationMethod const& integration_method,
    MeshLib::Element const& element, bool const is_axially_symmetric)
{
    auto const shape_matrices_at_element_center =
        NumLib::initShapeMatricesAtElementCenter<
            ShapeFunction, ShapeMatricesType, DisplacementDim>(
            element, is_axially_symmetric, integration_method);

    auto const N_0 = shape_matrices_at_element_center.N;
    auto const dNdx_0 = shape_matrices_at_element_center.dNdx;

    // For the 2D case the 33-component is needed (and the four entries
    // of the non-symmetric matrix); In 3d there are nine entries.
    GradientMatrixType G0(
        DisplacementDim * DisplacementDim + (DisplacementDim == 2 ? 1 : 0),
        ShapeFunction::NPOINTS * DisplacementDim);

    auto const x_coord =
        NumLib::interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
            element, N_0);
    Deformation::computeGMatrix<DisplacementDim, ShapeFunction::NPOINTS>(
        dNdx_0, G0, is_axially_symmetric, N_0, x_coord);

    GradientVectorType const grad_u = G0 * u;
    GradientVectorType const F0 =
        grad_u + MathLib::VectorizedTensor::identity<DisplacementDim>();

    if (compute_detF0_only)
    {
        return {MathLib::VectorizedTensor::determinant(F0), {}};
    }

    // Note: MathLib::VectorizedTensor::toTensor returns a (3 x 3) matrix by
    //       Eigen::Matrix3D. The size of the return matrix is wrong, which
    //       should be (DisplacementDim x DisplacementDim).
    // Besides, Eigen::Matrix member inverse does not work with the mapped
    // matrix, e.g. that by using topLeftCorner. Therefore a hard copy is
    // needed for the following inverse, i.e. declaring F0inv as type
    // Eigen::MatrixXd instead of type auto.
    //
    Eigen::MatrixXd const F0_matrix =
        MathLib::VectorizedTensor::toTensor<DisplacementDim>(F0)
            .template topLeftCorner<DisplacementDim, DisplacementDim>();
    auto const F0inv = F0_matrix.inverse();

    GlobalDimNodalMatrixType F0InvN =
        GlobalDimNodalMatrixType::Zero(DisplacementDim, ShapeFunction::NPOINTS);

    for (unsigned i = 0; i < ShapeFunction::NPOINTS; ++i)
    {
        auto const dNidx_0 = dNdx_0.col(i);

        for (int k = 0; k < DisplacementDim; k++)
        {
            auto const Finv_k_col = F0inv.col(k);
            F0InvN(k, i) = Finv_k_col.transpose() * dNidx_0;
        }
    }

    return {MathLib::VectorizedTensor::determinant(F0), F0InvN};
}

template <int DisplacementDim, int NPOINTS, typename GlobalDimNodalMatrixType,
          typename GradientVectorType, typename BMatrixType, typename DNDX_Type>
BMatrixType computeBFbarMatrix(
    double const alpha, DNDX_Type const& dNdx, GradientVectorType const& F,
    GlobalDimNodalMatrixType const& F0InvN,
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const& eps,
    bool const is_axially_symmetric)
{
    BMatrixType BFbar = BMatrixType::Zero(
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim),
        NPOINTS * DisplacementDim);

    // Note: MathLib::VectorizedTensor::toTensor returns a (3 x 3) matrix by
    //       Eigen::Matrix3D. The size of the return matrix is wrong, which
    //       should be (DisplacementDim x DisplacementDim).
    // Besides, Eigen::Matrix member inverse does not work with the mapped
    // matrix, e.g. that by using topLeftCorner. Therefore a hard copy is
    // needed for the following inverse, i.e. declaring F0inv as type
    // Eigen::MatrixXd instead of type auto.
    //
    Eigen::MatrixXd const F_matrix =
        MathLib::VectorizedTensor::toTensor<DisplacementDim>(F)
            .template topLeftCorner<DisplacementDim, DisplacementDim>();
    auto const Finv = F_matrix.inverse();

    double const fac = alpha * alpha / 3;
    for (int i = 0; i < NPOINTS; ++i)
    {
        auto const dNidx = dNdx.col(i);
        for (int k = 0; k < DisplacementDim; k++)
        {
            // K_{col k}
            //
            // (F^{-1})_{nk} in vector form
            auto const Finv_k_col = Finv.col(k);
            // (F^{-1})_{nk}\cdot \nabla N^{i}
            double const Finv_nk_grad_N_i =
                F0InvN(k, i) - Finv_k_col.transpose() * dNidx;

            auto BFbar_col = BFbar.col(i + k * NPOINTS);
            BFbar_col = 2.0 * eps;
            BFbar_col.template segment<3>(0) += Eigen::Vector3d::Constant(1.0);
            BFbar_col *= fac * Finv_nk_grad_N_i;

            if constexpr (DisplacementDim == 2)
            {
                BFbar_col(2) = 0.0;
            }
        }
    }

    (void)is_axially_symmetric;  // Not considered yet
    return BFbar;
}
}  // namespace NonLinearFbar
}  // namespace ProcessLib
