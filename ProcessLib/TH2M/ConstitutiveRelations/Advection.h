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
#include "ConstitutiveDensity.h"
#include "PermeabilityData.h"
#include "PhaseTransitionData.h"
#include "PureLiquidDensity.h"
#include "Saturation.h"
#include "Viscosity.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
struct AdvectionData
{
    GlobalDimMatrix<DisplacementDim> advection_C_G;
    GlobalDimMatrix<DisplacementDim> advection_C_L;
    GlobalDimMatrix<DisplacementDim> advection_W_G;
    GlobalDimMatrix<DisplacementDim> advection_W_L;
};

template <int DisplacementDim>
struct AdvectionDerivativeData
{
    GlobalDimMatrix<DisplacementDim> dadvection_C_dp_GR;
    GlobalDimMatrix<DisplacementDim> dadvection_C_dp_cap;
};

template <int DisplacementDim>
struct AdvectionModel
{
    void eval(ConstituentDensityData const& constituent_density_data,
              PermeabilityData<DisplacementDim> const& permeability_data,
              PureLiquidDensityData const& rho_W_LR,
              ViscosityData const& viscosity_data,
              AdvectionData<DisplacementDim>& advection_data) const;

    void dEval(
        ConstituentDensityData const& constituent_density_data,
        PermeabilityData<DisplacementDim> const& permeability_data,
        ViscosityData const& viscosity_data,
        SaturationDataDeriv const& dS_L_dp_cap,
        PhaseTransitionData const& phase_transition_data,
        AdvectionDerivativeData<DisplacementDim>& advection_d_data) const;
};

extern template struct AdvectionModel<2>;
extern template struct AdvectionModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
