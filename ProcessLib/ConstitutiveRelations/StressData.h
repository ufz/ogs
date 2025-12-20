// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
