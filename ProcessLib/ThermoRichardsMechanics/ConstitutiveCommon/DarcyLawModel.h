// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "CapillaryPressureData.h"
#include "DarcyLawData.h"
#include "LiquidDensityData.h"
#include "LiquidViscosityData.h"
#include "PermeabilityData.h"
#include "ProcessLib/ConstitutiveRelations/Base.h"
#include "ThermoOsmosisData.h"

namespace ProcessLib::ThermoRichardsMechanics
{

template <int DisplacementDim>
struct DarcyLawModel
{
    explicit DarcyLawModel(ProcessLib::ConstitutiveRelations::SpecificBodyForce<
                           DisplacementDim> const& specific_body_force)
        : specific_body_force_(specific_body_force)
    {
    }

    void eval(CapillaryPressureData<DisplacementDim> const& p_cap_data,
              LiquidDensityData const& rho_L_data,
              LiquidViscosityData const& mu_L_data,
              PermeabilityData<DisplacementDim> const& perm_data,
              ThermoOsmosisData<DisplacementDim> const& th_osmosis_data,
              DarcyLawData<DisplacementDim>& out) const;

    static DarcyLawModel create(
        ProcessLib::ConstitutiveRelations::SpecificBodyForce<
            DisplacementDim> const& specific_body_force)
    {
        return DarcyLawModel{specific_body_force};
    }

private:
    /// Gravity vector (specific body force).
    ProcessLib::ConstitutiveRelations::SpecificBodyForce<DisplacementDim> const
        specific_body_force_;
};

extern template struct DarcyLawModel<2>;
extern template struct DarcyLawModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
