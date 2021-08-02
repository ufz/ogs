/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <memory>
#include <vector>

namespace MaterialPropertyLib
{
class Property;
}

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
std::unique_ptr<Property> createTemperatureDependentDiffusion(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);
}  // namespace MaterialPropertyLib
