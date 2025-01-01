/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <memory>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace ParameterLib
{
struct ParameterBase;
}

namespace MaterialPropertyLib
{
class Property;
}

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createCubicLawPermeability(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);
}  // namespace MaterialPropertyLib
