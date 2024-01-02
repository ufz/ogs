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
struct LiquidViscosityData
{
    double viscosity;

    static auto reflect()
    {
        return ProcessLib::Reflection::reflectWithName(
            "viscosity", &LiquidViscosityData::viscosity);
    }
};

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
