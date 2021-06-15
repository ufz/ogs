/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include <memory>

namespace ProcessLib
{
namespace ThermoRichardsFlow
{
struct SimplifiedElasticityModel;
}
}

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
