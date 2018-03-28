/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "MaterialLib/SolidModels/MechanicsBase.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
template <typename ShapeMatrixType>
struct IntegrationPointDataSoil final
{
    explicit IntegrationPointDataSoil() {}

    typename ShapeMatrixType::NodalRowVectorType N;
    typename ShapeMatrixType::GlobalDimNodalMatrixType dNdx;
    double integration_weight;

    void pushBackState() {}

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
