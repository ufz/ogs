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
struct StrainData
{
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> eps;

    static auto reflect()
    {
        using Self = StrainData<DisplacementDim>;

        return ProcessLib::Reflection::reflectWithName("epsilon", &Self::eps);
    }
};

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
