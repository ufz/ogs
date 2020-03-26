/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 25, 2020, 5:00 PM
 */

#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialPropertyLib
{
class CapillaryPressureBrooksCorey;
}

namespace MaterialPropertyLib
{
std::unique_ptr<CapillaryPressureBrooksCorey>
createCapillaryPressureBrooksCorey(BaseLib::ConfigTree const& config);
}  // namespace MaterialPropertyLib
