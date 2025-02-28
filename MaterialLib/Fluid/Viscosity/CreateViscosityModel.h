/**
 * \file
 *  \brief A function for creating viscosity model
 *
 *  \copyright
 *   Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "BaseLib/ConfigTree-fwd.h"
#include "MaterialLib/Fluid/FluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
///  Create a viscosity model
///  \param config  ConfigTree object has a tag of `<viscosity>`
std::unique_ptr<FluidProperty> createViscosityModel(
    BaseLib::ConfigTree const& config);

}  // end namespace
}  // namespace MaterialLib
