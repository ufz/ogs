/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "Base.h"

namespace ProcessLib::ThermoRichardsMechanics
{

template <int DisplacementDim>
struct TotalStressData
{
    KelvinVector<DisplacementDim> sigma_total = KVnan<DisplacementDim>();

    static auto reflect()
    {
        using Self = TotalStressData<DisplacementDim>;

        return ProcessLib::Reflection::reflectWithName("sigma_total",
                                                       &Self::sigma_total);
    }
};
}  // namespace ProcessLib::ThermoRichardsMechanics
