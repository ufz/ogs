// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "ConstitutiveDensity.h"
#include "Enthalpy.h"
#include "FluidDensity.h"
#include "MassMoleFractions.h"
#include "PhaseTransitionData.h"
#include "PureLiquidDensity.h"
#include "VapourPartialPressure.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
struct PhaseTransitionModel
{
    explicit PhaseTransitionModel(
        std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const&
            media)
    {
        DBUG("Create phase transition models...");

        // check for minimum requirement definitions in media object
        std::array const required_gas_properties = {
            MaterialPropertyLib::density};
        std::array const required_liquid_properties = {
            MaterialPropertyLib::density};

        for (auto const& m : media)
        {
            checkRequiredProperties(m.second->phase("Gas"),
                                    required_gas_properties);
            checkRequiredProperties(m.second->phase("AqueousLiquid"),
                                    required_liquid_properties);
        }
    }

    virtual ~PhaseTransitionModel() = default;

    virtual void eval(SpaceTimeData const& x_t, MediaData const& media_data,
                      GasPressureData const& p_GR,
                      CapillaryPressureData const& p_cap,
                      TemperatureData const& T_data,
                      PureLiquidDensityData const& rho_W_LR,
                      FluidEnthalpyData& fluid_enthalpy_data,
                      MassMoleFractionsData& mass_mole_fractions_data,
                      FluidDensityData& fluid_density_data,
                      VapourPartialPressureData& vapour_pressure_data,
                      ConstituentDensityData& constituent_density_data,
                      PhaseTransitionData& cv) const = 0;
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
