// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "LiquidDensityData.h"
#include "LiquidViscosityData.h"
#include "MediaData.h"
#include "TemperatureData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
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
