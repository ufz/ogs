// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"

namespace ProcessLib::ThermoRichardsMechanics
{
struct LiquidDensityData
{
    double rho_LR;
    double drho_LR_dp;
    double drho_LR_dT;

    static auto reflect()
    {
        return ProcessLib::Reflection::reflectWithName(
            "liquid_density", &LiquidDensityData::rho_LR);
    }
};

template <int DisplacementDim>
struct LiquidDensityModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              CapillaryPressureData<DisplacementDim> const& p_cap_data,
              TemperatureData<DisplacementDim> const& T_data,
              LiquidDensityData& out) const;
};

extern template struct LiquidDensityModel<2>;
extern template struct LiquidDensityModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
