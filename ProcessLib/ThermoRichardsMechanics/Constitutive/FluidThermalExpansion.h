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

#include "LiquidDensity.h"
#include "MathLib/KelvinVector.h"
#include "Porosity.h"
#include "SolidThermalExpansion.h"

namespace ProcessLib::ThermoRichardsMechanics
{
struct FluidThermalExpansionData
{
    double eff_thermal_expansion;
};

template <int DisplacementDim>
struct FluidThermalExpansionModel
{
    void eval(
        SpaceTimeData const& x_t, MediaData const& media_data,
        SolidThermalExpansionData<DisplacementDim> const& s_therm_exp_data,
        PorosityData const& poro_data, LiquidDensityData const& rho_L_data,
        BiotData const& biot_data, FluidThermalExpansionData& out) const;
};

extern template struct FluidThermalExpansionModel<2>;
extern template struct FluidThermalExpansionModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics