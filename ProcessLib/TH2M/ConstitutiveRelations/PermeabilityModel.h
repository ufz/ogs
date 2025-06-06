/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Base.h"
#include "EquivalentPlasticStrainData.h"
#include "PermeabilityData.h"
#include "ProcessLib/ConstitutiveRelations/StrainData.h"
#include "Saturation.h"
#include "TotalStress.h"
#include "TransportPorosity.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
struct PermeabilityModel
{
    void eval(
        SpaceTimeData const& x_t, MediaData const& media_data,
        SaturationData const& S_L_data, GasPressureData const& p_g,
        CapillaryPressureData const& p_cap, TemperatureData const& T_data,
        TransportPorosityData const& transport_poro_data,
        TotalStressData<DisplacementDim> const& total_stress_data,
        MechanicalStrainData<DisplacementDim> const& mechanical_strain_data,
        StrainData<DisplacementDim> const& eps_data,
        EquivalentPlasticStrainData const& equivalent_plastic_strain,
        PermeabilityData<DisplacementDim>& out) const;
};

extern template struct PermeabilityModel<2>;
extern template struct PermeabilityModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
