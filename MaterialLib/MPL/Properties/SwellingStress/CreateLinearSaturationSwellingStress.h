/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on August 14, 2020, 8:47 AM
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
std::unique_ptr<Property> createLinearSaturationSwellingStress(
    BaseLib::ConfigTree const& config);
}  // namespace MaterialPropertyLib
