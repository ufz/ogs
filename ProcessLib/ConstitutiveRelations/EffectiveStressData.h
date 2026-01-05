// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
