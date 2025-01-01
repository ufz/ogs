/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
