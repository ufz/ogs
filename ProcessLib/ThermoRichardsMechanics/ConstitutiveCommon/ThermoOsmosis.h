/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
