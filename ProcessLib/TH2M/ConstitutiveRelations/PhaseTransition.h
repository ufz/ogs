// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <map>

#include "MaterialLib/MPL/Medium.h"
#include "PhaseTransitionModel.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
struct PhaseTransition : PhaseTransitionModel
{
    explicit PhaseTransition(
        std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const&
            media);

    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              GasPressureData const& p_GR, CapillaryPressureData const& p_cap,
              TemperatureData const& T_data,
              PureLiquidDensityData const& rho_W_LR,
              FluidEnthalpyData& fluid_enthalpy_data,
              MassMoleFractionsData& mass_mole_fractions_data,
              FluidDensityData& fluid_density_data,
              VapourPartialPressureData& vapour_pressure_data,
              ConstituentDensityData& constituent_density_data,
              PhaseTransitionData& cv) const override;

private:
    int const n_components_gas_;
    int const gas_phase_vapour_component_index_;
    int const gas_phase_dry_air_component_index_;
    int const liquid_phase_solute_component_index_;
    int const liquid_phase_solvent_component_index_;
};

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
