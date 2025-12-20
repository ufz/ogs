// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "BaseLib/StrongType.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{

using SolidHeatCapacityData =
    BaseLib::StrongType<double, struct SolidHeatCapacityTag>;

struct SolidHeatCapacityModel
{
    void eval(SpaceTimeData const& x_t,
              MediaData const& media_data,
              TemperatureData const& T_data,
              SolidHeatCapacityData& solid_heat_capacity) const;
};

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
