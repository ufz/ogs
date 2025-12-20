// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>

namespace ProcessLib
{
namespace ThermoRichardsFlow
{
struct SimplifiedElasticityModel;
}
}  // namespace ProcessLib

namespace BaseLib
{
class ConfigTree;
}

namespace ProcessLib
{
namespace ThermoRichardsFlow
{
std::unique_ptr<SimplifiedElasticityModel> createElasticityModel(
    BaseLib::ConfigTree const& config);

}  // namespace ThermoRichardsFlow
}  // namespace ProcessLib
