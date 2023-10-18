/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Base.h"

namespace ProcessLib::SmallDeformation
{

template <int DisplacementDim>
struct StressData
{
    // Stress is stateful for some constitutive settings and therefore must be
    // initialized to something valid, e.g., zero.
    // TODO find a better solution for that.
    KelvinVector<DisplacementDim> sigma = KVzero<DisplacementDim>();

    static auto reflect()
    {
        using Self = StressData<DisplacementDim>;

        return ProcessLib::Reflection::reflectWithName("sigma", &Self::sigma);
    }
};
}  // namespace ProcessLib::SmallDeformation
