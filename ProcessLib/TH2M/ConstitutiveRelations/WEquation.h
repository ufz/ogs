/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Advection.h"
#include "Base.h"
#include "Biot.h"
#include "ConstitutiveDensity.h"
#include "FluidDensity.h"
#include "PermeabilityData.h"
#include "PhaseTransitionData.h"
#include "Porosity.h"
#include "Saturation.h"
#include "SolidCompressibility.h"
#include "Viscosity.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
struct FW1Data
{
    GlobalDimMatrix<DisplacementDim> A;
};

template <int DisplacementDim>
struct FW1Model
{
    void eval(AdvectionData<DisplacementDim> const& advection_data,
              FluidDensityData const& fluid_density_data,
              FW1Data<DisplacementDim>& fW_1) const;
};

extern template struct FW1Model<2>;
extern template struct FW1Model<3>;

struct FW2Data
{
    double a = nan;
};

struct FW2DerivativeData
{
    double dp_GR = nan;
    double dp_cap = nan;
    double dT = nan;
};

struct FW2Model
{
    void eval(BiotData const biot_data,
              CapillaryPressureData const pCap,
              ConstituentDensityData const& constituent_density_data,
              PorosityData const& porosity_data,
              PureLiquidDensityData const& rho_W_LR,
              SaturationData const& S_L_data,
              SolidCompressibilityData const beta_p_SR,
              FW2Data& fW_2) const;

    void dEval(BiotData const& biot_data,
               CapillaryPressureData const pCap,
               ConstituentDensityData const& constituent_density_data,
               PhaseTransitionData const& phase_transition_data,
               PorosityData const& porosity_data,
               PorosityDerivativeData const& porosity_d_data,
               PureLiquidDensityData const& rho_W_LR,
               SaturationData const& S_L_data,
               SaturationDataDeriv const& dS_L_dp_cap,
               SolidCompressibilityData const& beta_p_SR,
               FW2DerivativeData& dfW_2) const;
};

struct FW3aData
{
    double a = nan;
};

struct FW3aDerivativeData
{
    double dp_GR = nan;
    double dp_cap = nan;
    double dT = nan;
};

struct FW3aModel
{
    void eval(
        double const dt,
        ConstituentDensityData const& constituent_density_data,
        PrevState<ConstituentDensityData> const& constituent_density_data_prev,
        PrevState<PureLiquidDensityData> const& rho_W_LR_prev,
        PureLiquidDensityData const& rho_W_LR,
        SaturationData const& S_L_data,
        FW3aData& fW_3a) const;

    void dEval(
        double const dt,
        ConstituentDensityData const& constituent_density_data,
        PhaseTransitionData const& phase_transition_data,
        PrevState<ConstituentDensityData> const& constituent_density_data_prev,
        PrevState<PureLiquidDensityData> const& rho_W_LR_prev,
        PureLiquidDensityData const& rho_W_LR,
        SaturationData const& S_L_data,
        SaturationDataDeriv const& dS_L_dp_cap,
        FW3aDerivativeData& dfW_3a) const;
};

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
