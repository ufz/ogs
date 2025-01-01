/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <map>

#include "MaterialLib/MPL/Medium.h"
#include "PhaseTransitionModel.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
struct NoPhaseTransition : PhaseTransitionModel
{
    explicit NoPhaseTransition(
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
};

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
