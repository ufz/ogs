/**
 * \file
 * \author Norbert Grunwald
 * \date   Jul 07 2020
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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

namespace MaterialPropertyLib
{
class AverageMolarMass;
}

namespace MaterialPropertyLib
{
std::unique_ptr<AverageMolarMass> createAverageMolarMass(
    BaseLib::ConfigTree const& config);
}  // namespace MaterialPropertyLib
