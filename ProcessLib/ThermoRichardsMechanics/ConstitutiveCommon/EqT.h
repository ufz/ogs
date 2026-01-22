// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "EqTData.h"
#include "TRMHeatStorageAndFluxData.h"
#include "TRMVaporDiffusionData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
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
