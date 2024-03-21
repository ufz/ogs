/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Base.h"
#include "PhaseTransitionData.h"
#include "PureLiquidDensity.h"
#include "Viscosity.h"

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
            MaterialPropertyLib::viscosity, MaterialPropertyLib::density};
        std::array const required_liquid_properties = {
            MaterialPropertyLib::viscosity, MaterialPropertyLib::density};

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
                      ViscosityData& viscosity_data,
                      PhaseTransitionData& cv) const = 0;
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
