/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "NumLib/Fem/AverageGradShapeFunction.h"
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
            FInvN[k * NPOINTS + i] = Finv.col(k).dot(dNidx);
        }
    }

    return FInvN;
}

template <int DisplacementDim, typename GradientVectorType,
          typename VectorTypeForFbar, typename NodalVectorType,
          typename ShapeFunction, typename ShapeMatricesType, typename IpData>
std::tuple<double, VectorTypeForFbar> computeFBarInitialVariablesAverage(
    std::vector<IpData, Eigen::aligned_allocator<IpData>> const& ip_data,
    bool const compute_detF0_only, Eigen::VectorXd const& u,
    NumLib::GenericIntegrationMethod const& integration_method,
    MeshLib::Element const& element, bool const is_axially_symmetric)
{
    unsigned const n_integration_points =
        integration_method.getNumberOfPoints();
    // Compute the element volume
    double volume = 0.0;
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& w = ip_data[ip].integration_weight;
        volume += w;
    }

    VectorTypeForFbar averaged_grad_N =
        VectorTypeForFbar::Zero(DisplacementDim * ShapeFunction::NPOINTS);
    NodalVectorType averaged_N_div_r;
    if (is_axially_symmetric)
    {
        averaged_N_div_r = NodalVectorType::Zero(ShapeFunction::NPOINTS);
    }

    GradientVectorType F0 =
        MathLib::VectorizedTensor::identity<DisplacementDim>();

    for (unsigned i = 0; i < ShapeFunction::NPOINTS; i++)
    {
        Eigen::Vector3d const bar_gradN =
            NumLib::averageGradShapeFunction<DisplacementDim, ShapeFunction,
                                             ShapeMatricesType, IpData>(
                i, element, integration_method, ip_data, is_axially_symmetric) /
            volume;
        averaged_grad_N.template segment<DisplacementDim>(i * DisplacementDim) =
            bar_gradN.template segment<DisplacementDim>(0);

        for (int k = 0; k < DisplacementDim; k++)
        {
            F0.template segment<DisplacementDim>(k * DisplacementDim) +=
                u[k * ShapeFunction::NPOINTS + i] *
                bar_gradN.template segment<DisplacementDim>(0);
        }
        if (is_axially_symmetric)
        {
            averaged_N_div_r[i] = bar_gradN[2];
            F0[4] += bar_gradN[2] * u[i];
        }
    }

    if (compute_detF0_only)
    {
        return {MathLib::VectorizedTensor::determinant(F0), {}};
    }

    VectorTypeForFbar F0InvN =
        VectorTypeForFbar::Zero(DisplacementDim * ShapeFunction::NPOINTS);

    Eigen::MatrixXd const F_matrix =
        MathLib::VectorizedTensor::toTensor<DisplacementDim>(F0)
            .template topLeftCorner<DisplacementDim, DisplacementDim>();
    auto const Finv = F_matrix.inverse();

    for (unsigned i = 0; i < ShapeFunction::NPOINTS; ++i)
    {
        auto const dNidx = averaged_grad_N.template segment<DisplacementDim>(
            i * DisplacementDim);
        for (int k = 0; k < DisplacementDim; k++)
        {
            F0InvN[k * ShapeFunction::NPOINTS + i] = Finv.col(k).dot(dNidx);
        }
        // TODO: if(is_axially_symmetric)
    }

    return {MathLib::VectorizedTensor::determinant(F0), F0InvN};
}

template <int DisplacementDim, typename GradientVectorType,
          typename GradientMatrixType, typename VectorTypeForFbar,
          typename ShapeFunction, typename ShapeMatricesType>
std::tuple<double, VectorTypeForFbar> computeFBarInitialVariablesElementCenter(
    bool const compute_detF0_only, Eigen::VectorXd const& u,
    MeshLib::Element const& element, bool const is_axially_symmetric)
{
    auto const shape_matrices_at_element_center =
        NumLib::initShapeMatricesAtElementCenter<
            ShapeFunction, ShapeMatricesType, DisplacementDim>(
            element, is_axially_symmetric);

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
