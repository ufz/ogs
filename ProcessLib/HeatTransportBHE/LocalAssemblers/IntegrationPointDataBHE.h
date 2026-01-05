// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace ProcessLib
{
namespace HeatTransportBHE
{
template <typename ShapeMatrixType>
struct IntegrationPointDataBHE final
{
    typename ShapeMatrixType::NodalRowVectorType const N;
    typename ShapeMatrixType::GlobalDimNodalMatrixType const dNdx;
    double const integration_weight;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
