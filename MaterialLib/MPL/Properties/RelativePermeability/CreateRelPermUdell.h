/**
 * \file
 * \author Norbert Grunwald
 * \date   01.12.2020
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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
class RelPermUdell;
}

namespace MaterialPropertyLib
{
std::unique_ptr<RelPermUdell> createRelPermUdell(
    BaseLib::ConfigTree const& config);
}  // namespace MaterialPropertyLib
