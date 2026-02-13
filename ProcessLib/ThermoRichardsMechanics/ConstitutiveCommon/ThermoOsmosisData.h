// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct ThermoOsmosisData
{
    GlobalDimMatrix<DisplacementDim> K_pT_Laplace;
    GlobalDimMatrix<DisplacementDim> K_Tp_Laplace;
    GlobalDimVector<DisplacementDim> seepage_velocity_contribution;
};
// Explicit instantiation declarations to avoid multiple-definition issues.
extern template struct ThermoOsmosisData<2>;
extern template struct ThermoOsmosisData<3>;

}  // namespace ProcessLib::ThermoRichardsMechanics
