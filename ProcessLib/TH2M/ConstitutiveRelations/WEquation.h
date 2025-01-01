/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
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

template <int DisplacementDim>
struct FW4LWpGData
{
    GlobalDimMatrix<DisplacementDim> L;
};

template <int DisplacementDim>
struct FW4LWpGDerivativeData
{
    GlobalDimMatrix<DisplacementDim> dp_GR;
    GlobalDimMatrix<DisplacementDim> dp_cap;
    GlobalDimMatrix<DisplacementDim> dT;
};

template <int DisplacementDim>
struct FW4LWpGModel
{
    void eval(AdvectionData<DisplacementDim> const& advection_data,
              FluidDensityData const& fluid_density_data,
              PhaseTransitionData const& phase_transition_data,
              PorosityData const& porosity_data,
              SaturationData const& S_L_data,
              FW4LWpGData<DisplacementDim>& fW_4_LWpG) const;

    void dEval(ConstituentDensityData const& constituent_density_data,
               PermeabilityData<DisplacementDim> const& permeability_data,
               PhaseTransitionData const& phase_transition_data,
               PureLiquidDensityData const& rho_W_LR,
               SaturationDataDeriv const& dS_L_dp_cap,
               ViscosityData const& viscosity_data,
               FW4LWpGDerivativeData<DisplacementDim>& dfW_4_LWpG) const;
};

extern template struct FW4LWpGModel<2>;
extern template struct FW4LWpGModel<3>;

template <int DisplacementDim>
struct FW4LWpCData
{
    GlobalDimMatrix<DisplacementDim> L;
};

template <int DisplacementDim>
struct FW4LWpCDerivativeData
{
    GlobalDimMatrix<DisplacementDim> dp_GR;
    GlobalDimMatrix<DisplacementDim> dp_cap;
    GlobalDimMatrix<DisplacementDim> dT;
};

template <int DisplacementDim>
struct FW4LWpCModel
{
    void eval(AdvectionData<DisplacementDim> const& advection_data,
              FluidDensityData const& fluid_density_data,
              PhaseTransitionData const& phase_transition_data,
              PorosityData const& porosity_data,
              SaturationData const& S_L_data,
              FW4LWpCData<DisplacementDim>& fW_4_LWpC) const;

    void dEval(AdvectionData<DisplacementDim> const& advection_data,
               FluidDensityData const& fluid_density_data,
               PermeabilityData<DisplacementDim> const& permeability_data,
               PhaseTransitionData const& phase_transition_data,
               PorosityData const& porosity_data,
               PureLiquidDensityData const& rho_W_LR,
               SaturationData const& S_L_data,
               SaturationDataDeriv const& dS_L_dp_cap,
               ViscosityData const& viscosity_data,
               FW4LWpCDerivativeData<DisplacementDim>& dfW_4_LWpC) const;
};

extern template struct FW4LWpCModel<2>;
extern template struct FW4LWpCModel<3>;

template <int DisplacementDim>
struct FW4LWTData
{
    GlobalDimMatrix<DisplacementDim> L;
};

template <int DisplacementDim>
struct FW4LWTModel
{
    void eval(FluidDensityData const& fluid_density_data,
              PhaseTransitionData const& phase_transition_data,
              PorosityData const& porosity_data,
              SaturationData const& S_L_data,
              FW4LWTData<DisplacementDim>& fW_4_LWT) const;
};

struct FW4MWpGData
{
    double m = nan;
};

struct FW4MWpGModel
{
    void eval(BiotData const& biot_data,
              ConstituentDensityData const& constituent_density_data,
              PorosityData const& porosity_data,
              PureLiquidDensityData const& rho_W_LR,
              SaturationData const& S_L_data,
              SolidCompressibilityData const& beta_p_SR,
              FW4MWpGData& fW_4_MWpG) const;
};

struct FW4MWpCData
{
    double m = nan;
    double ml = nan;
};

struct FW4MWpCModel
{
    void eval(BiotData const& biot_data,
              CapillaryPressureData const pCap,
              ConstituentDensityData const& constituent_density_data,
              PorosityData const& porosity_data,
              PrevState<SaturationData> const& S_L_data_prev,
              PureLiquidDensityData const& rho_W_LR,
              SaturationData const& S_L_data,
              SolidCompressibilityData const& beta_p_SR,
              FW4MWpCData& fW_4_MWpC) const;
};

struct FW4MWTData
{
    double m = nan;
};

template <int DisplacementDim>
struct FW4MWTModel
{
    void eval(
        BiotData const& biot_data,
        ConstituentDensityData const& constituent_density_data,
        PorosityData const& porosity_data,
        PureLiquidDensityData const& rho_W_LR,
        SaturationData const& S_L_data,
        SolidThermalExpansionData<DisplacementDim> const& s_therm_exp_data,
        FW4MWTData& fW_4_MWT) const;
};

extern template struct FW4MWTModel<2>;
extern template struct FW4MWTModel<3>;

struct FW4MWuData
{
    double m = nan;
};

struct FW4MWuModel
{
    void eval(BiotData const& biot_data,
              ConstituentDensityData const& constituent_density_data,
              PureLiquidDensityData const& rho_W_LR,
              SaturationData const& S_L_data,
              FW4MWuData& fW_4_MWu) const;
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
