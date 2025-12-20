// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "Saturation.h"

namespace ProcessLib::ThermoRichardsMechanics
{
struct BishopsData
{
    double chi_S_L = nan;
    double dchi_dS_L = nan;
};

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
