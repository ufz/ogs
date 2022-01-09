/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on January 26, 2021, 11:11 AM
 */

#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialPropertyLib
{
class RelPermUdellNonwettingPhase;
}

namespace MaterialPropertyLib
{
std::unique_ptr<RelPermUdellNonwettingPhase> createRelPermUdellNonwettingPhase(
    BaseLib::ConfigTree const& config);
}  // namespace MaterialPropertyLib
