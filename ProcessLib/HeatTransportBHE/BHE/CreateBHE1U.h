/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BHE_1U.h"
#include "MaterialLib/Fluid/FluidProperty.h"
#include "ProcessLib/Process.h"
namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE  // namespace of borehole heat exchanger
{
BHE::BHE_1U* CreateBHE1U(
    BaseLib::ConfigTree const& bhe_conf,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> const&
        bhe_refrigerant_density,
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> const&
        bhe_refrigerant_viscosity,
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> const&
        bhe_refrigerant_heat_capacity,
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> const&
        bhe_regrigerant_heat_conductivity);
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
