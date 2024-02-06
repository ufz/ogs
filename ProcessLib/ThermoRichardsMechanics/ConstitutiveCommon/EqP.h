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

#include "Biot.h"
#include "FluidThermalExpansion.h"
#include "LiquidDensity.h"
#include "LiquidViscosity.h"
#include "PermeabilityData.h"
#include "Saturation.h"
#include "TRMStorage.h"
#include "TRMVaporDiffusion.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct EqPData
{
    GlobalDimVector<DisplacementDim> J_pp_dNT_V_N = DVnan<DisplacementDim>();
    double J_pp_X_BTI2NT_u_dot_N = nan;

    GlobalDimMatrix<DisplacementDim> K_pp_Laplace = DMnan<DisplacementDim>();

    double M_pT_X_NTN = nan;
    double M_pu_X_BTI2N = nan;

    GlobalDimVector<DisplacementDim> rhs_p_dNT_V = DVnan<DisplacementDim>();

    double storage_p_a_p_X_NTN = nan;
};

template <int DisplacementDim>
struct EqPModel
{
    explicit EqPModel(
        Eigen::Vector<double, DisplacementDim> const& specific_body_force)
        : b_(specific_body_force)
    {
    }

    void eval(CapillaryPressureData<DisplacementDim> const& p_cap_data,
              TemperatureData<DisplacementDim> const& T_data,
              SaturationData const& S_L_data,
              SaturationDataDeriv const& dS_L_data,
              BiotData const& biot_data,
              LiquidDensityData const& rho_L_data,
              LiquidViscosityData const& mu_L_data,
              PermeabilityData<DisplacementDim> const& perm_data,
              FluidThermalExpansionData const& f_therm_exp_data,
              TRMVaporDiffusionData<DisplacementDim> const& vap_data,
              TRMStorageData const& storage_data,
              EqPData<DisplacementDim>& out) const;

private:
    /// Gravity vector (specific body force).
    Eigen::Vector<double, DisplacementDim> const b_;
};

extern template struct EqPModel<2>;
extern template struct EqPModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
