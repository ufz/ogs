/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialPropertyLib
{
class SaturationExponential;
}

namespace MaterialPropertyLib
{
std::unique_ptr<SaturationExponential> createSaturationExponential(
    BaseLib::ConfigTree const& config);
}  // namespace MaterialPropertyLib
