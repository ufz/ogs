/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Base.h"
#include "BaseLib/StrongType.h"
#include "ProcessLib/Reflection/ReflectionData.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
struct SaturationData
{
    double S_L = nan;

    static auto reflect()
    {
        return ProcessLib::Reflection::reflectWithName("saturation",
                                                       &SaturationData::S_L);
    }
};

using SaturationDataDeriv =
    BaseLib::StrongType<double, struct SaturationDataDerivTag>;

struct SaturationModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              CapillaryPressureData const& p_cap, SaturationData& S_L_data,
              SaturationDataDeriv& dS_L_data) const;
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
