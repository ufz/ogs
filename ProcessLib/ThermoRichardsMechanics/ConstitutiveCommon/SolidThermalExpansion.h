// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "MediaData.h"
#include "SolidThermalExpansionData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct SolidThermalExpansionModel
{
    void eval(SpaceTimeData const& x_t,
              MediaData const& media_data,
              SolidThermalExpansionData<DisplacementDim>& out) const;
};

extern template struct SolidThermalExpansionModel<2>;
extern template struct SolidThermalExpansionModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
