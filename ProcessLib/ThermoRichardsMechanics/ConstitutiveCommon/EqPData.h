// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct EqPData
{
    GlobalDimVector<DisplacementDim> J_pp_dNT_V_N = DVnan<DisplacementDim>();
    double J_pp_X_BTI2NT_u_dot_N = nan;

    GlobalDimMatrix<DisplacementDim> K_pp_Laplace = DMnan<DisplacementDim>();

    double M_pT_X_NTN = nan;
    double M_pu_X_BTI2N = nan;

    GlobalDimVector<DisplacementDim> rhs_p_dNT_V = DVnan<DisplacementDim>();

    double storage_p_a_p_X_NTN = nan;
};
// Explicit instantiation declarations to avoid multiple-definition issues.
extern template struct EqPData<2>;
extern template struct EqPData<3>;

}  // namespace ProcessLib::ThermoRichardsMechanics
