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
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
