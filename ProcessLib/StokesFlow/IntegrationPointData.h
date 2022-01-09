/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace ProcessLib::StokesFlow
{
template <typename ShapeMatrixTypeLiquidVelocity,
          typename ShapeMatrixTypePressure,
          int GlobalDim,
          int NPoints>
struct IntegrationPointData final
{
    using NodalRowVectorTypeLiquidVelocity =
        typename ShapeMatrixTypeLiquidVelocity::NodalRowVectorType;
    using GlobalDimNodalMatrixTypeLiquidVelocity =
        typename ShapeMatrixTypeLiquidVelocity::GlobalDimNodalMatrixType;

    using NodalRowVectorTypePressure =
        typename ShapeMatrixTypePressure::NodalRowVectorType;
    using GlobalDimNodalMatrixTypePressure =
        typename ShapeMatrixTypePressure::GlobalDimNodalMatrixType;

    IntegrationPointData(NodalRowVectorTypeLiquidVelocity const& N_v_,
                         GlobalDimNodalMatrixTypeLiquidVelocity const& dNdx_v_,
                         NodalRowVectorTypePressure const& N_p_,
                         GlobalDimNodalMatrixTypePressure const& dNdx_p_,
                         double const& integration_weight_)
        : N_v(N_v_),
          dNdx_v(dNdx_v_),
          N_p(N_p_),
          dNdx_p(dNdx_p_),
          integration_weight(integration_weight_)
    {
    }

    NodalRowVectorTypeLiquidVelocity const N_v;
    GlobalDimNodalMatrixTypeLiquidVelocity const dNdx_v;

    typename ShapeMatrixTypeLiquidVelocity::
        template MatrixType<GlobalDim, NPoints * GlobalDim>
            N_v_op;
    typename ShapeMatrixTypeLiquidVelocity::
        template MatrixType<GlobalDim * GlobalDim, NPoints * GlobalDim>
            dNdx_v_op;

    NodalRowVectorTypePressure const N_p;
    GlobalDimNodalMatrixTypePressure const dNdx_p;

    double const integration_weight;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};
}  // namespace ProcessLib::StokesFlow
