// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace ProcessLib::LargeDeformation
{

template <typename BMatricesType, typename ShapeMatricesType,
          int DisplacementDim>
struct IntegrationPointData final
{
    double integration_weight;
    typename ShapeMatricesType::NodalRowVectorType N_u;
    typename ShapeMatricesType::GlobalDimNodalMatrixType dNdx_u;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ProcessLib::LargeDeformation
