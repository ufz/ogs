// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "BaseLib/StrongType.h"
#include "LiquidDensity.h"

namespace ProcessLib::ThermoRichardsMechanics
{
using LiquidViscosityData =
    BaseLib::StrongType<double, struct LiquidViscosityDataTag>;

constexpr std::string_view ioName(struct LiquidViscosityDataTag*)
{
    return "viscosity";
}

template <int DisplacementDim>
struct LiquidViscosityModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
            LiquidDensityData const& rho_L_data,
            TemperatureData<DisplacementDim> const& T_data,
              LiquidViscosityData& out) const;
};

extern template struct LiquidViscosityModel<2>;
extern template struct LiquidViscosityModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
