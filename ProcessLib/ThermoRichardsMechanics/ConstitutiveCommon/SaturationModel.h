// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "CapillaryPressureData.h"
#include "MediaData.h"
#include "SaturationData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct SaturationModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              CapillaryPressureData<DisplacementDim> const& p_cap_data,
              SaturationData& S_L_data, SaturationDataDeriv& dS_L_data) const;
};

extern template struct SaturationModel<2>;
extern template struct SaturationModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
