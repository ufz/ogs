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

}  // namespace ProcessLib::ThermoRichardsMechanics
