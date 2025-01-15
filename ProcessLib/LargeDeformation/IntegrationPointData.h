/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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
