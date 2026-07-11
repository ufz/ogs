// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct TRMVaporDiffusionData
{
    double heat_capacity_vapor = 0;
    GlobalDimVector<DisplacementDim> vapor_flux;
    double storage_coefficient_by_water_vapor = 0;

    double J_pT_X_dNTdN = 0;
    double K_pp_X_dNTdN = 0;
    double K_TT_X_dNTdN = 0;
    double K_Tp_X_dNTdN = 0;
    double M_Tp_X_NTN = 0;
    double M_TT_X_NTN = 0;
    double M_pT_X_NTN = 0;

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
// Explicit instantiation declarations to avoid multiple-definition issues.
extern template struct TRMVaporDiffusionData<2>;
extern template struct TRMVaporDiffusionData<3>;

}  // namespace ProcessLib::ThermoRichardsMechanics
