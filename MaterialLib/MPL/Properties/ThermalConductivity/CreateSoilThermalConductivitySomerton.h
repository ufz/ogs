/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on February 17, 2021, 3:47 PM
 */

#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialPropertyLib
{
class Property;

std::unique_ptr<Property> createSoilThermalConductivitySomerton(
    BaseLib::ConfigTree const& config);
}  // namespace MaterialPropertyLib
