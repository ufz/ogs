// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"

namespace ProcessLib::ThermoRichardsMechanics
{

template <int DisplacementDim>
struct TotalStressData
{
    // Total stress is stateful for some constitutive settings and therefore
    // must be initialized to something valid, e.g., zero.
    // TODO find a better solution for that.
    KelvinVector<DisplacementDim> sigma_total = KV::KVzero<DisplacementDim>();

    static auto reflect()
    {
        using Self = TotalStressData<DisplacementDim>;

        return ProcessLib::Reflection::reflectWithName("sigma_total",
                                                       &Self::sigma_total);
    }
};
}  // namespace ProcessLib::ThermoRichardsMechanics
