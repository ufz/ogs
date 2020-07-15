/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 */

#pragma once

#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

namespace ProcessLib
{

template <typename NodalRowVectorType>
struct SourceTermIntegrationPointData final
{
    SourceTermIntegrationPointData(NodalRowVectorType N_,
                                   double const& integration_weight_)
        : N(std::move(N_)), integration_weight(integration_weight_)
    {
    }

    NodalRowVectorType const N;
    double const integration_weight;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ProcessLib
