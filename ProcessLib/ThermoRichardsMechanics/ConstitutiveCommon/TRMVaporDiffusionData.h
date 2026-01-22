// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct TRMVaporDiffusionData
{
    double heat_capacity_vapor;
    GlobalDimVector<DisplacementDim> vapor_flux;
    double storage_coefficient_by_water_vapor;

    double J_pT_X_dNTdN;
    double K_pp_X_dNTdN;
    double K_TT_X_dNTdN;
    double K_Tp_X_dNTdN;
    double M_Tp_X_NTN;
    double M_TT_X_NTN;
    double M_pT_X_NTN;

    void setZero()
    {
        heat_capacity_vapor = 0;
        vapor_flux = GlobalDimVector<DisplacementDim>::Zero(DisplacementDim);
        storage_coefficient_by_water_vapor = 0;

        J_pT_X_dNTdN = 0;
        K_pp_X_dNTdN = 0;
        K_TT_X_dNTdN = 0;
        K_Tp_X_dNTdN = 0;
        M_Tp_X_NTN = 0;
        M_TT_X_NTN = 0;
        M_pT_X_NTN = 0;
    }
};

}  // namespace ProcessLib::ThermoRichardsMechanics
