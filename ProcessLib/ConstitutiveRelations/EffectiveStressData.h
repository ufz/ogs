/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
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
struct EffectiveStressData
{
    // TODO it seems fragile that some data have to be initialized that way.
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> sigma_eff =
        MathLib::KelvinVector::KVzero<DisplacementDim>();

    static auto reflect()
    {
        using Self = EffectiveStressData<DisplacementDim>;

        return ProcessLib::Reflection::reflectWithName("sigma",
                                                       &Self::sigma_eff);
    }
};

}  // namespace ProcessLib::ConstitutiveRelations
