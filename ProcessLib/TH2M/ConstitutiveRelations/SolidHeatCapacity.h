/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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
