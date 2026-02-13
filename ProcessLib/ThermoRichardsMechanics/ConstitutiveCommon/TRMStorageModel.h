// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Biot.h"
#include "CapillaryPressureData.h"
#include "LiquidDensityData.h"
#include "PorosityData.h"
#include "ProcessLib/ConstitutiveRelations/Base.h"
#include "SaturationData.h"
#include "SolidCompressibilityData.h"
#include "TRMStorageData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
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
