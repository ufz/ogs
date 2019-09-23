/**
 *  \brief A function for creating a thermal conductivity model for fluid
 *
 *  \copyright
 *   Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *   \file
 */

#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialLib
{
namespace Fluid
{
class FluidProperty;

/**
 *  Create a thermal conductivity model
 *  \param config  ConfigTree object has a tag of `<thermal_conductivity>`
 */
std::unique_ptr<FluidProperty> createFluidThermalConductivityModel(
    BaseLib::ConfigTree const& config);

}  // namespace Fluid
}  // namespace MaterialLib
