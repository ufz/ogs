/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on April 21, 2020, 11:06 AM
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
std::unique_ptr<Property> createCapillaryPressureRegularizedVanGenuchten(
    BaseLib::ConfigTree const& config);
}  // namespace MaterialPropertyLib
