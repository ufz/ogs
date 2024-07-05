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

#include <boost/algorithm/string/predicate.hpp>
#include <cmath>
#include <string_view>
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

enum class BarDetFType
{
    ELEMENT_CENTER_VALUE,
    ELEMENT_AVERAGE,
    NONE
};

BarDetFType convertStringToDetFBarType(
    std::string_view const bar_det_f_type_name);

template <int DisplacementDim, int NPOINTS, typename VectorTypeForFbar,
          typename GradientVectorType, typename DNDX_Type>
VectorTypeForFbar computeQVector(DNDX_Type const& dNdx,
                                 GradientVectorType const& F,
                                 bool const /*is_axially_symmetric*/)
{
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

    VectorTypeForFbar FInvN =
        VectorTypeForFbar::Zero(DisplacementDim * NPOINTS);
    for (int i = 0; i < NPOINTS; ++i)
    {
        auto const dNidx = dNdx.col(i);
        for (int k = 0; k < DisplacementDim; k++)
        {
            auto const Finv_k_col = Finv.col(k);

            FInvN[k * NPOINTS + i] = Finv_k_col.dot(dNidx);
        }
    }

    return FInvN;
}

template <int DisplacementDim, typename GradientVectorType,
          typename GradientMatrixType, typename VectorTypeForFbar,
          typename ShapeFunction, typename ShapeMatricesType>
std::tuple<double, VectorTypeForFbar> computeFBarInitialVariables(
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

    VectorTypeForFbar F0InvN =
        computeQVector<DisplacementDim, ShapeFunction::NPOINTS,
                       VectorTypeForFbar, GradientVectorType,
                       typename ShapeMatricesType::GlobalDimNodalMatrixType>(
            dNdx_0, F0, is_axially_symmetric);

    return {MathLib::VectorizedTensor::determinant(F0), F0InvN};
}

template <int DisplacementDim>
MathLib::KelvinVector::KelvinVectorType<DisplacementDim> identityForF(
    bool const is_axially_symmetric)
{
    if (DisplacementDim == 3 || is_axially_symmetric)
    {
        return MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::kelvin_vector_dimensions(
                DisplacementDim)>::identity2;
    }

    auto identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    identity2[2] = 0.0;

    return identity2;
}

template <int DisplacementDim>
MathLib::KelvinVector::KelvinVectorType<DisplacementDim> compute2EPlusI(
    double const alpha_p2,
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const& eps_bar,
    bool const is_axially_symmetric)
{
    return (2.0 * eps_bar +
            identityForF<DisplacementDim>(is_axially_symmetric)) /
           alpha_p2;
}

}  // namespace NonLinearFbar
}  // namespace ProcessLib
