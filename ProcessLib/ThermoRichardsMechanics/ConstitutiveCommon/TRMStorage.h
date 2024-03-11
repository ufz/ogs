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

#include "Biot.h"
#include "LiquidDensity.h"
#include "Porosity.h"
#include "Saturation.h"

namespace ProcessLib::ThermoRichardsMechanics
{
struct TRMStorageData
{
    double storage_p_a_p;
    double storage_p_a_S_X_NTN;
    double J_pp_X_NTN;
    double storage_p_a_S_Jpp_X_NTN;
};

template <int DisplacementDim>
struct TRMStorageModel
{
    void eval(SpaceTimeData const& x_t, BiotData const& biot_data,
              PorosityData const& poro_data,
              LiquidDensityData const& rho_L_data,
              SaturationData const& S_L_data,
              SaturationDataDeriv const& dS_L_data,
              PrevState<SaturationData> const& S_L_prev_data,
              CapillaryPressureData<DisplacementDim> const& p_cap_data,
              SolidCompressibilityData const& solid_compressibility_data,
              TRMStorageData& out) const;
};

extern template struct TRMStorageModel<2>;
extern template struct TRMStorageModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
