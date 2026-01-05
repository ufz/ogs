// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "TRMHeatStorageAndFlux.h"
#include "TRMVaporDiffusion.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct EqTData
{
    GlobalDimVector<DisplacementDim> K_TT_NT_V_dN = DVnan<DisplacementDim>();
    double M_TT_X_NTN = nan;
};

template <int DisplacementDim>
struct EqTModel
{
    void eval(TRMHeatStorageAndFluxData<DisplacementDim> const& heat_data,
              TRMVaporDiffusionData<DisplacementDim> const& vap_data,
              EqTData<DisplacementDim>& out) const;
};

extern template struct EqTModel<2>;
extern template struct EqTModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
