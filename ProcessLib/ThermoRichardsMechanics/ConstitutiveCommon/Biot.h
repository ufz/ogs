// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "BaseLib/StrongType.h"

namespace ProcessLib::ThermoRichardsMechanics
{
using BiotData = BaseLib::StrongType<double, struct BiotTag>;

struct BiotModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              BiotData& out) const;
};
}  // namespace ProcessLib::ThermoRichardsMechanics
