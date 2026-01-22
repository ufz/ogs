// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct TRMHeatStorageAndFluxData
{
    double M_TT_X_NTN;
    GlobalDimMatrix<DisplacementDim> K_TT_Laplace;
    GlobalDimVector<DisplacementDim> K_Tp_NT_V_dN;
    double K_Tp_X_NTN;
    GlobalDimVector<DisplacementDim>
        advective_heat_flux_contribution_to_K_liquid;
};

}  // namespace ProcessLib::ThermoRichardsMechanics
