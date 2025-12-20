// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "LiquidDensity.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct ThermoOsmosisData
{
    GlobalDimMatrix<DisplacementDim> K_pT_Laplace;
    GlobalDimMatrix<DisplacementDim> K_Tp_Laplace;
    GlobalDimVector<DisplacementDim> seepage_velocity_contribution;
};

template <int DisplacementDim>
struct ThermoOsmosisModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              TemperatureData<DisplacementDim> const& T_data,
              LiquidDensityData const& rho_L_data,
              ThermoOsmosisData<DisplacementDim>& out) const;
};

extern template struct ThermoOsmosisModel<2>;
extern template struct ThermoOsmosisModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
