// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct GravityData
{
    GlobalDimVector<DisplacementDim> volumetric_body_force;
    GlobalDimVector<DisplacementDim> J_up_HT_V_N;
};
// Explicit instantiation declarations to avoid multiple-definition issues.
extern template struct GravityData<2>;
extern template struct GravityData<3>;

}  // namespace ProcessLib::ThermoRichardsMechanics
