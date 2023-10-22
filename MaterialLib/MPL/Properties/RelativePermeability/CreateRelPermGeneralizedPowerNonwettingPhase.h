/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on January 26, 2024, 11:11 AM
 */

#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialPropertyLib
{
class RelPermGeneralizedPowerNonwettingPhase;
}

namespace MaterialPropertyLib
{
std::unique_ptr<RelPermGeneralizedPowerNonwettingPhase>
createRelPermGeneralizedPowerNonwettingPhase(BaseLib::ConfigTree const& config);
}  // namespace MaterialPropertyLib
