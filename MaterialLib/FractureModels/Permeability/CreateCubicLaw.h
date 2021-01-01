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

namespace BaseLib
{
class ConfigTree;
}
namespace MaterialLib::Fracture::Permeability
{
class CubicLaw;
class Permeability;
}  // namespace MaterialLib::Fracture::Permeability

namespace MaterialLib::Fracture::Permeability
{
std::unique_ptr<Permeability> createCubicLaw(BaseLib::ConfigTree const& config);
}  // namespace MaterialLib::Fracture::Permeability
