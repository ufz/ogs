// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MathLib/KelvinVector.h"
#include "ProcessLib/Reflection/ReflectionData.h"

namespace ProcessLib::ConstitutiveRelations
{
template <int DisplacementDim>
struct MechanicalStrainData
{
    // TODO it seems fragile that some data have to be initialized that way.
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> eps_m =
        MathLib::KelvinVector::KVzero<DisplacementDim>();

    static auto reflect()
    {
        using Self = MechanicalStrainData<DisplacementDim>;

        return ProcessLib::Reflection::reflectWithName("eps_m", &Self::eps_m);
    }
};

}  // namespace ProcessLib::ConstitutiveRelations
