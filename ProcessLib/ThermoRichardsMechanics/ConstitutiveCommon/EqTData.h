// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct EqTData
{
    GlobalDimVector<DisplacementDim> K_TT_NT_V_dN = DVnan<DisplacementDim>();
    double M_TT_X_NTN = nan;
};
// Explicit instantiation declarations to avoid multiple-definition issues.
extern template struct EqTData<2>;
extern template struct EqTData<3>;

}  // namespace ProcessLib::ThermoRichardsMechanics
