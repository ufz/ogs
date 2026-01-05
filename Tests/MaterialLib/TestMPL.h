// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
