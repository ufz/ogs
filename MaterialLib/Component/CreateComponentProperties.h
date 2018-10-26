/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   CreateComponentProperties.h
 */

#pragma once

#include <memory>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace ProcessLib
{
struct ParameterBase;
}

namespace MaterialLib
{
namespace Component
{
class ComponentProperties;

std::vector<ComponentProperties> createComponentProperties(
    BaseLib::ConfigTree const& configs,
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters);
}  // namespace Component
}  // namespace MaterialLib
