// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "CapillaryPressureData.h"
#include "EquivalentPlasticStrainData.h"
#include "MediaData.h"
#include "PermeabilityData.h"
#include "ProcessLib/ConstitutiveRelations/StrainData.h"
#include "SaturationData.h"
#include "TemperatureData.h"
#include "TotalStressData.h"
#include "TransportPorosityData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct PermeabilityModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              SaturationData const& S_L_data,
              CapillaryPressureData<DisplacementDim> const& p_cap_data,
              TemperatureData<DisplacementDim> const& T_data,
              TransportPorosityData const& transport_poro_data,
              TotalStressData<DisplacementDim> const& total_stress_data,
              StrainData<DisplacementDim> const& eps_data,
              EquivalentPlasticStrainData const& equiv_plast_strain_data,
              PermeabilityData<DisplacementDim>& out) const;
};

extern template struct PermeabilityModel<2>;
extern template struct PermeabilityModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
