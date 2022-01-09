/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <memory>

#include "Permeability.h"

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialLib::Fracture::Permeability
{
std::unique_ptr<Permeability> createPermeabilityModel(
    BaseLib::ConfigTree const& config);
}  // namespace MaterialLib::Fracture::Permeability
