/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "Porosity.h"

namespace ProcessLib::ThermoRichardsMechanics
{
struct LiquidViscosityData
{
    double viscosity;
};

template <int DisplacementDim>
struct LiquidViscosityModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              TemperatureData<DisplacementDim> const& T_data,
              LiquidViscosityData& out) const;
};

extern template struct LiquidViscosityModel<2>;
extern template struct LiquidViscosityModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics