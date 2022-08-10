/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "LiquidViscosity.h"
#include "Porosity.h"
#include "SolidMechanics.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct PermeabilityData
{
    double k_rel;
    double dk_rel_dS_L;
    GlobalDimMatrix<DisplacementDim> Ki_over_mu;
};

template <int DisplacementDim>
struct PermeabilityModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              SolidCompressibilityData const& solid_compressibility_data,
              SaturationData const& S_L_data, BishopsData const& bishops_data,
              BishopsData const& bishops_data_prev,
              CapillaryPressureData<DisplacementDim> const& p_cap_data,
              TemperatureData<DisplacementDim> const& T_data,
              PorosityData const& poro_data,
              LiquidViscosityData const& mu_L_data,
              // TODO evaluate transport porosity evolution separately
              TransportPorosityData& transport_poro_data,
              TransportPorosityData const& transport_poro_data_prev,
              SolidMechanicsDataStateless<DisplacementDim> const& s_mech_data,
              StrainData<DisplacementDim> const& eps_data,
              StrainData<DisplacementDim> const& eps_prev_data,
              PermeabilityData<DisplacementDim>& out) const;
};

extern template struct PermeabilityModel<2>;
extern template struct PermeabilityModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
