/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on January 25, 2021, 11:49 AM
 */

#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialPropertyLib
{
class RelPermBrooksCoreyNonwettingPhase;
}

namespace MaterialPropertyLib
{
std::unique_ptr<RelPermBrooksCoreyNonwettingPhase>
createRelPermBrooksCoreyNonwettingPhase(BaseLib::ConfigTree const& config);
}  // namespace MaterialPropertyLib
