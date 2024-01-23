/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "MathLib/KelvinVector.h"
#include "ProcessLib/Reflection/ReflectionData.h"

namespace ProcessLib::ConstitutiveRelations
{
template <int DisplacementDim>
struct StressData
{
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> sigma;

    static auto reflect()
    {
        using Self = StressData<DisplacementDim>;

        return ProcessLib::Reflection::reflectWithName("sigma", &Self::sigma);
    }
};
}  // namespace ProcessLib::ConstitutiveRelations
