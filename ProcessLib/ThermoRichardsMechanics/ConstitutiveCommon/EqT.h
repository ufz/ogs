/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
