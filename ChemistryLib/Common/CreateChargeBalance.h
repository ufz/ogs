// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace BaseLib
{
class ConfigTree;
}

namespace ChemistryLib
{
enum class ChargeBalance;

ChargeBalance createChargeBalance(BaseLib::ConfigTree const& config);
}  // namespace ChemistryLib
