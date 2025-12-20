// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MathLib/KelvinVector.h"
#include "ProcessLib/Reflection/ReflectionData.h"

namespace ProcessLib::ConstitutiveRelations
{
template <int DisplacementDim>
struct StrainData
{
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> eps;

    static auto reflect()
    {
        using Self = StrainData<DisplacementDim>;

        return ProcessLib::Reflection::reflectWithName("epsilon", &Self::eps);
    }
};

}  // namespace ProcessLib::ConstitutiveRelations
