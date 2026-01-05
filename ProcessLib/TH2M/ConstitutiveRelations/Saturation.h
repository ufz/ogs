// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
              CapillaryPressureData const& p_cap,
              SaturationData& S_L_data) const;

    void dEval(SpaceTimeData const& x_t, MediaData const& media_data,
               CapillaryPressureData const& p_cap,
               SaturationDataDeriv& dS_L_data) const;
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
