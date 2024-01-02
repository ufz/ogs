/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Base.h"
#include "BaseLib/StrongType.h"

namespace ProcessLib::SmallDeformation
{
using SolidDensity = BaseLib::StrongType<double, struct SolidDensityTag>;

struct SolidDensityModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              Temperature const& temperature, SolidDensity& out) const;
};
}  // namespace ProcessLib::SmallDeformation
