// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Biot.h"
#include "CapillaryPressureData.h"
#include "EqPData.h"
#include "FluidThermalExpansionData.h"
#include "LiquidDensityData.h"
#include "LiquidViscosityData.h"
#include "PermeabilityData.h"
#include "ProcessLib/ConstitutiveRelations/Base.h"
#include "SaturationData.h"
#include "TRMStorageData.h"
#include "TRMVaporDiffusionData.h"
#include "TemperatureData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct EqPModel
{
    explicit EqPModel(ProcessLib::ConstitutiveRelations::SpecificBodyForce<
                      DisplacementDim> const& specific_body_force)
        : specific_body_force_(specific_body_force)
    {
    }

    void eval(CapillaryPressureData<DisplacementDim> const& p_cap_data,
              TemperatureData<DisplacementDim> const& T_data,
              SaturationData const& S_L_data,
              SaturationDataDeriv const& dS_L_data,
              BiotData const& biot_data,
              LiquidDensityData const& rho_L_data,
              LiquidViscosityData const& mu_L_data,
              PermeabilityData<DisplacementDim> const& perm_data,
              FluidThermalExpansionData const& f_therm_exp_data,
              TRMVaporDiffusionData<DisplacementDim> const& vap_data,
              TRMStorageData const& storage_data,
              EqPData<DisplacementDim>& out) const;

    static EqPModel create(ProcessLib::ConstitutiveRelations::SpecificBodyForce<
                           DisplacementDim> const& specific_body_force)
    {
        return EqPModel{specific_body_force};
    }

private:
    /// Gravity vector (specific body force).
    ProcessLib::ConstitutiveRelations::SpecificBodyForce<DisplacementDim> const
        specific_body_force_;
};

extern template struct EqPModel<2>;
extern template struct EqPModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
