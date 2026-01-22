// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "BishopsData.h"
#include "MediaData.h"
#include "SaturationData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
struct BishopsModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              SaturationData const& S_L_data, BishopsData& out) const;
};

struct BishopsPrevModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              PrevState<SaturationData> const& S_L_data,
              PrevState<BishopsData>& out) const;
};
}  // namespace ProcessLib::ThermoRichardsMechanics
