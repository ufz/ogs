/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   CreateFluidProperties.h
 *
 * Created on December 13, 2016, 3:32 PM
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
class FluidProperties;
/// Create an instance of class FluidProperties
/// \param config  ConfigTree object has tags of `<fluid>`
std::unique_ptr<FluidProperties> createFluidProperties(
    BaseLib::ConfigTree const& config);
}  // end namespace
}  // end namespace
