// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
struct FC1Data
{
    GlobalDimMatrix<DisplacementDim> A;
};

template <int DisplacementDim>
struct FC1Model
{
    void eval(AdvectionData<DisplacementDim> const& advection_data,
              FluidDensityData const& fluid_density_data,
              FC1Data<DisplacementDim>& fC_1) const;
};

extern template struct FC1Model<2>;
extern template struct FC1Model<3>;

struct FC2aData
{
    double a = nan;
};

struct FC2aDerivativeData
{
    double dp_GR = nan;
    double dp_cap = nan;
    double dT = nan;
};

struct FC2aModel
{
    void eval(BiotData const biot_data,
              CapillaryPressureData const pCap,
              ConstituentDensityData const& constituent_density_data,
              PorosityData const& porosity_data,
              SaturationData const& S_L_data,
              SolidCompressibilityData const beta_p_SR,
              FC2aData& fC_2a) const;

    void dEval(BiotData const& biot_data,
               CapillaryPressureData const pCap,
               ConstituentDensityData const& constituent_density_data,
               PhaseTransitionData const& phase_transition_data,
               PorosityData const& porosity_data,
               PorosityDerivativeData const& porosity_d_data,
               SaturationData const& S_L_data,
               SaturationDataDeriv const& dS_L_dp_cap,
               SolidCompressibilityData const& beta_p_SR,
               FC2aDerivativeData& dfC_2a) const;
};

struct FC3aData
{
    double a = nan;
};

struct FC3aDerivativeData
{
    double dp_GR = nan;
    double dp_cap = nan;
    double dT = nan;
};

struct FC3aModel
{
    void eval(
        double const dt,
        ConstituentDensityData const& constituent_density_data,
        PrevState<ConstituentDensityData> const& constituent_density_data_prev,
        SaturationData const& S_L_data,
        FC3aData& fC_3a) const;

    void dEval(
        double const dt,
        ConstituentDensityData const& constituent_density_data,
        PrevState<ConstituentDensityData> const& constituent_density_data_prev,
        PhaseTransitionData const& phase_transition_data,
        SaturationData const& S_L_data,
        SaturationDataDeriv const& dS_L_dp_cap,
        FC3aDerivativeData& dfC_3a) const;
};

template <int DisplacementDim>
struct FC4LCpGData
{
    GlobalDimMatrix<DisplacementDim> L;
};

template <int DisplacementDim>
struct FC4LCpGDerivativeData
{
    GlobalDimMatrix<DisplacementDim> dp_GR;
    GlobalDimMatrix<DisplacementDim> dp_cap;
    GlobalDimMatrix<DisplacementDim> dT;
};

template <int DisplacementDim>
struct FC4LCpGModel
{
    void eval(AdvectionData<DisplacementDim> const& advection_data,
              FluidDensityData const& fluid_density_data,
              PhaseTransitionData const& phase_transition_data,
              PorosityData const& porosity_data,
              SaturationData const& S_L_data,
              FC4LCpGData<DisplacementDim>& fC_4_LCpG) const;

    void dEval(PermeabilityData<DisplacementDim> const& permeability_data,
               ViscosityData const& viscosity_data,
               PhaseTransitionData const& phase_transition_data,
               AdvectionDerivativeData<DisplacementDim> const& advection_d_data,
               FC4LCpGDerivativeData<DisplacementDim>& dfC_4_LCpG) const;
};

extern template struct FC4LCpGModel<2>;
extern template struct FC4LCpGModel<3>;

template <int DisplacementDim>
struct FC4LCpCData
{
    GlobalDimMatrix<DisplacementDim> L;
};

template <int DisplacementDim>
struct FC4LCpCDerivativeData
{
    GlobalDimMatrix<DisplacementDim> dp_GR;
    GlobalDimMatrix<DisplacementDim> dp_cap;
    GlobalDimMatrix<DisplacementDim> dT;
};

template <int DisplacementDim>
struct FC4LCpCModel
{
    void eval(AdvectionData<DisplacementDim> const& advection_data,
              FluidDensityData const& fluid_density_data,
              PhaseTransitionData const& phase_transition_data,
              PorosityData const& porosity_data,
              SaturationData const& S_L_data,
              FC4LCpCData<DisplacementDim>& fC_4_LCpC) const;

    void dEval(ConstituentDensityData const& constituent_density_data,
               PermeabilityData<DisplacementDim> const& permeability_data,
               PhaseTransitionData const& phase_transition_data,
               SaturationDataDeriv const& dS_L_dp_cap,
               ViscosityData const& viscosity_data,
               FC4LCpCDerivativeData<DisplacementDim>& dfC_4_LCpC) const;
};

extern template struct FC4LCpCModel<2>;
extern template struct FC4LCpCModel<3>;

template <int DisplacementDim>
struct FC4LCTData
{
    GlobalDimMatrix<DisplacementDim> L;
};

template <int DisplacementDim>
struct FC4LCTModel
{
    void eval(FluidDensityData const& fluid_density_data,
              PhaseTransitionData const& phase_transition_data,
              PorosityData const& porosity_data,
              SaturationData const& S_L_data,
              FC4LCTData<DisplacementDim>& fC_4_LCT) const;
};

struct FC4MCpGData
{
    double m = nan;
};

struct FC4MCpGDerivativeData
{
    double dp_GR = nan;
    double dT = nan;
};

struct FC4MCpGModel
{
    void eval(BiotData const& biot_data,
              ConstituentDensityData const& constituent_density_data,
              PorosityData const& porosity_data,
              SaturationData const& S_L_data,
              SolidCompressibilityData const& beta_p_SR,
              FC4MCpGData& fC_4_MCpG) const;

    void dEval(BiotData const& biot_data,
               ConstituentDensityData const& constituent_density_data,
               PhaseTransitionData const& phase_transition_data,
               PorosityData const& porosity_data,
               PorosityDerivativeData const& porosity_d_data,
               SaturationData const& S_L_data,
               SolidCompressibilityData const& beta_p_SR,
               FC4MCpGDerivativeData& dfC_4_MCpG) const;
};

struct FC4MCpCData
{
    double m = nan;
    double ml = nan;
};

struct FC4MCpCModel
{
    void eval(BiotData const& biot_data,
              CapillaryPressureData const pCap,
              ConstituentDensityData const& constituent_density_data,
              PorosityData const& porosity_data,
              PrevState<SaturationData> const& S_L_data_prev,
              SaturationData const& S_L_data,
              SolidCompressibilityData const& beta_p_SR,
              FC4MCpCData& fC_4_MCpC) const;
};

struct FC4MCTData
{
    double m = nan;
};

struct FC4MCTDerivativeData
{
    double dT = nan;
};

template <int DisplacementDim>
struct FC4MCTModel
{
    void eval(
        BiotData const& biot_data,
        ConstituentDensityData const& constituent_density_data,
        PorosityData const& porosity_data,
        SaturationData const& S_L_data,
        SolidThermalExpansionData<DisplacementDim> const& s_therm_exp_data,
        FC4MCTData& fC_4_MCT) const;

    void dEval(
        BiotData const& biot_data,
        ConstituentDensityData const& constituent_density_data,
        PhaseTransitionData const& phase_transition_data,
        PorosityData const& porosity_data,
        PorosityDerivativeData const& porosity_d_data,
        SaturationData const& S_L_data,
        SolidThermalExpansionData<DisplacementDim> const& s_therm_exp_data,
        FC4MCTDerivativeData& dfC_4_MCT) const;
};

struct FC4MCuData
{
    double m = nan;
};

struct FC4MCuDerivativeData
{
    double dT = nan;
};

struct FC4MCuModel
{
    void eval(BiotData const& biot_data,
              ConstituentDensityData const& constituent_density_data,
              SaturationData const& S_L_data,
              FC4MCuData& fC_4_MCu) const;

    void dEval(BiotData const& biot_data,
               PhaseTransitionData const& phase_transition_data,
               SaturationData const& S_L_data,
               FC4MCuDerivativeData& dfC_4_MCu) const;
};

extern template struct FC4MCTModel<2>;
extern template struct FC4MCTModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
