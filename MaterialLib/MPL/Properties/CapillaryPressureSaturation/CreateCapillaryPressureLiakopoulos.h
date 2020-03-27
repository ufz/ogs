/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 27, 2020, 12:27 PM
 */

#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialPropertyLib
{
class CapillaryPressureLiakopoulos;
}

namespace MaterialPropertyLib
{
std::unique_ptr<CapillaryPressureLiakopoulos>
createCapillaryPressureLiakopoulos(BaseLib::ConfigTree const& config);
}  // namespace MaterialPropertyLib
