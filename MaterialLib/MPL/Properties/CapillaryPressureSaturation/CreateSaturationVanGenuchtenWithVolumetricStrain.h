// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialPropertyLib
{
class SaturationVanGenuchtenWithVolumetricStrain;
}

namespace MaterialPropertyLib
{
std::unique_ptr<SaturationVanGenuchtenWithVolumetricStrain>
createSaturationVanGenuchtenWithVolumetricStrain(
    BaseLib::ConfigTree const& config);
}  // namespace MaterialPropertyLib
