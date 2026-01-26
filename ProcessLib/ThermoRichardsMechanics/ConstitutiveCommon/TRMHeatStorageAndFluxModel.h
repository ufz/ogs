// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "DarcyLawData.h"
#include "LiquidDensityData.h"
#include "LiquidViscosityData.h"
#include "MediaData.h"
#include "PermeabilityData.h"
#include "PorosityData.h"
#include "SaturationData.h"
#include "SolidDensityData.h"
#include "TRMHeatStorageAndFluxData.h"
#include "TemperatureData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct TRMHeatStorageAndFluxModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              LiquidDensityData const& rho_L_data,
              SolidDensityData const& rho_S_data,
              SaturationData const& S_L_data,
              SaturationDataDeriv const& dS_L_data,
              PorosityData const& poro_data,
              LiquidViscosityData const& mu_L_data,
              PermeabilityData<DisplacementDim> const& perm,
              TemperatureData<DisplacementDim> const& T_data,
              DarcyLawData<DisplacementDim> const& darcy_data,
              TRMHeatStorageAndFluxData<DisplacementDim>& out) const;
};

extern template struct TRMHeatStorageAndFluxModel<2>;
extern template struct TRMHeatStorageAndFluxModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
