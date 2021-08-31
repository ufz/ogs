/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "MathLib/KelvinVector.h"

namespace ProcessLib::TH2M
{
/// Variables needed only for the assembly process. The values are not preserved
/// throughout the iterations contrary to the variables in IntegrationPointData.
template <int DisplacementDim>
struct ConstitutiveVariables
{
    using KelvinMatrixType =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    KelvinMatrixType C;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ProcessLib::TH2M
