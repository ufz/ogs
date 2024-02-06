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

#include "ProcessLib/ConstitutiveRelations/StrainData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/EquivalentPlasticStrainData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/PermeabilityData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TotalStressData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TransportPorosity.h"

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
