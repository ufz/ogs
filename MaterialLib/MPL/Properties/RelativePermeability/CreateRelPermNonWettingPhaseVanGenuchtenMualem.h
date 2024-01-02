/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on April 14, 2020, 3:50 PM
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
std::unique_ptr<Property> createRelPermNonWettingPhaseVanGenuchtenMualem(
    BaseLib::ConfigTree const& config);
}  // namespace MaterialPropertyLib
