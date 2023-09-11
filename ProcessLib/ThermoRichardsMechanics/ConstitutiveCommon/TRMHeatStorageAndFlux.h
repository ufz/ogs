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

#include "DarcyLaw.h"
#include "LiquidDensity.h"
#include "PermeabilityData.h"
#include "Porosity.h"
#include "SolidDensity.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct TRMHeatStorageAndFluxData
{
    double M_TT_X_NTN;
    GlobalDimMatrix<DisplacementDim> K_TT_Laplace;
    GlobalDimVector<DisplacementDim> K_Tp_NT_V_dN;
    double K_Tp_X_NTN;
    GlobalDimVector<DisplacementDim>
        advective_heat_flux_contribution_to_K_liquid;
};

template <int DisplacementDim>
struct TRMHeatStorageAndFluxModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              LiquidDensityData const& rho_L_data,
              SolidDensityData const& rho_S_data,
              SaturationData const& S_L_data,
              SaturationDataDeriv const& dS_L_data,
              PorosityData const& poro_data,
              LiquidViscosityData const& mu_L_data,
              PermeabilityData<DisplacementDim> const& perm,
              TemperatureData<DisplacementDim> const& T_data,
              DarcyLawData<DisplacementDim> const& darcy_data,
              TRMHeatStorageAndFluxData<DisplacementDim>& out) const;
};

extern template struct TRMHeatStorageAndFluxModel<2>;
extern template struct TRMHeatStorageAndFluxModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
