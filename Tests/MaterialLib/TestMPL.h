/**
 * \file
 * \author Norbert Grunwald
 * \date   Oct 22, 2018
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <functional>
#include <memory>

#include "MaterialLib/MPL/Medium.h"

namespace MPL = MaterialPropertyLib;

namespace ParameterLib
{
struct CoordinateSystem;
}

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialPropertyLib
{
class Property;
}

namespace Tests
{
std::unique_ptr<MPL::Medium> createTestMaterial(
    std::string const& xml,
    int const geometry_dimension = 1,
    ParameterLib::CoordinateSystem const* const local_coordinate_system =
        nullptr);

std::unique_ptr<MaterialPropertyLib::Property> createTestProperty(
    const char xml[],
    std::function<std::unique_ptr<MaterialPropertyLib::Property>(
        BaseLib::ConfigTree const& config)>
        createProperty);

std::string makeConstantPropertyElement(std::string const name,
                                        double const value);

}  // namespace Tests
