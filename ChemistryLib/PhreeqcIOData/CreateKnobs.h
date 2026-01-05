// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace BaseLib
{
class ConfigTree;
}

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct Knobs;

Knobs createKnobs(BaseLib::ConfigTree const& config);
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
