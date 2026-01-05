// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "Porosity.h"
#include "ProcessLib/ConstitutiveRelations/EffectiveStressData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Bishops.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveStress_StrainTemperature/SolidMechanics.h"

namespace ProcessLib::ThermoRichardsMechanics
{
struct SolidDensityData
{
    double rho_SR;
    double dry_density_solid;

    static auto reflect()
    {
        return ProcessLib::Reflection::reflectWithName(
            "dry_density_solid", &SolidDensityData::dry_density_solid);
    }
};

template <int DisplacementDim>
struct SolidDensityModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              PorosityData const& poro_data,
              TemperatureData<DisplacementDim> const& T_data,
              ProcessLib::ConstitutiveRelations::EffectiveStressData<
                  DisplacementDim> const& sigma_eff_data,
              CapillaryPressureData<DisplacementDim> const& p_cap_data,
              BishopsData const& bishops_data, SolidDensityData& out) const;
};

extern template struct SolidDensityModel<2>;
extern template struct SolidDensityModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
