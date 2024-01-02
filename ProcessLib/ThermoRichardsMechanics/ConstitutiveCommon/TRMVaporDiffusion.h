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

#include "LiquidDensity.h"
#include "Porosity.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct TRMVaporDiffusionData
{
    double heat_capacity_vapor;
    GlobalDimVector<DisplacementDim> vapor_flux;
    double storage_coefficient_by_water_vapor;

    double J_pT_X_dNTdN;
    double K_pp_X_dNTdN;
    double K_TT_X_dNTdN;
    double K_Tp_X_dNTdN;
    double M_Tp_X_NTN;
    double M_TT_X_NTN;
    double M_pT_X_NTN;

    void setZero();
};

template <int DisplacementDim>
struct TRMVaporDiffusionModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              LiquidDensityData const& rho_L_data,
              SaturationData const& S_L_data,
              SaturationDataDeriv const& dS_L_data,
              PorosityData const& poro_data,
              CapillaryPressureData<DisplacementDim> const& p_cap_data,
              TemperatureData<DisplacementDim> const& T_data,
              TRMVaporDiffusionData<DisplacementDim>& out) const;
};

extern template struct TRMVaporDiffusionData<2>;
extern template struct TRMVaporDiffusionData<3>;
extern template struct TRMVaporDiffusionModel<2>;
extern template struct TRMVaporDiffusionModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
