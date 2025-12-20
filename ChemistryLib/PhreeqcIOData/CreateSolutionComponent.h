// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct Component;

std::vector<Component> createSolutionComponents(
    BaseLib::ConfigTree const& config);
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
