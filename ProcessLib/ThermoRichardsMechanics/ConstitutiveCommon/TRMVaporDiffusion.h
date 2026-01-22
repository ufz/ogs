// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "CapillaryPressureData.h"
#include "LiquidDensityData.h"
#include "MediaData.h"
#include "PorosityData.h"
#include "SaturationData.h"
#include "TRMVaporDiffusionData.h"
#include "TemperatureData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct TRMVaporDiffusionModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              LiquidDensityData const& rho_L_data,
              SaturationData const& S_L_data,
              SaturationDataDeriv const& dS_L_data,
              PorosityData const& poro_data,
              CapillaryPressureData<DisplacementDim> const& p_cap_data,
              TemperatureData<DisplacementDim> const& T_data,
              TRMVaporDiffusionData<DisplacementDim>& out) const;
};

extern template struct TRMVaporDiffusionModel<2>;
extern template struct TRMVaporDiffusionModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
