// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
