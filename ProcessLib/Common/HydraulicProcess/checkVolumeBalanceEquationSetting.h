// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string_view>

namespace MaterialPropertyLib
{
class MaterialSpatialDistributionMap;
}

namespace ProcessLib::Common::HydraulicProcess
{

void checkVolumeBalanceEquationSetting(
    MaterialPropertyLib::MaterialSpatialDistributionMap const& media_map);

}  // namespace ProcessLib::Common::HydraulicProcess
