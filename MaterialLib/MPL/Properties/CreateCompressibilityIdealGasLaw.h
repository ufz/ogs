/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on January 28, 2020, 16:05 PM
 */

#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialPropertyLib
{
class CompressibilityIdealGasLaw;
}

namespace MaterialPropertyLib
{
std::unique_ptr<CompressibilityIdealGasLaw> createCompressibilityIdealGasLaw(
    BaseLib::ConfigTree const& config);
}  // namespace MaterialPropertyLib
