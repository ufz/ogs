// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "CapillaryPressureData.h"
#include "LiquidDensityData.h"
#include "MediaData.h"
#include "TemperatureData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct LiquidDensityModel
{
    void eval(SpaceTimeData const& x_t,
              MediaData const& media_data,
              CapillaryPressureData<DisplacementDim> const& p_cap_data,
              TemperatureData<DisplacementDim> const& T_data,
              LiquidDensityData& out) const;
};

extern template struct LiquidDensityModel<2>;
extern template struct LiquidDensityModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
