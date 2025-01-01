/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CEquation.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
void FC1Model<DisplacementDim>::eval(
    AdvectionData<DisplacementDim> const& advection_data,
    FluidDensityData const& fluid_density_data,
    FC1Data<DisplacementDim>& fC_1) const
{
    fC_1.A = advection_data.advection_C_G * fluid_density_data.rho_GR +
             advection_data.advection_C_L * fluid_density_data.rho_LR;
}

template struct FC1Model<2>;
template struct FC1Model<3>;

void FC2aModel::eval(BiotData const biot_data,
                     CapillaryPressureData const pCap,
                     ConstituentDensityData const& constituent_density_data,
                     PorosityData const& porosity_data,
                     SaturationData const& S_L_data,
                     SolidCompressibilityData const beta_p_SR,
                     FC2aData& fC_2a) const
{
    auto const S_L = S_L_data.S_L;
    auto const S_G = 1. - S_L;
    double const rho_C_FR = S_G * constituent_density_data.rho_C_GR +
                            S_L * constituent_density_data.rho_C_LR;
    fC_2a.a =
        porosity_data.phi * (constituent_density_data.rho_C_LR -
                             constituent_density_data.rho_C_GR) -
        rho_C_FR * pCap() * (biot_data() - porosity_data.phi) * beta_p_SR();
}

void FC2aModel::dEval(BiotData const& biot_data,
                      CapillaryPressureData const pCap,
                      ConstituentDensityData const& constituent_density_data,
                      PhaseTransitionData const& phase_transition_data,
                      PorosityData const& porosity_data,
                      PorosityDerivativeData const& porosity_d_data,
                      SaturationData const& S_L_data,
                      SaturationDataDeriv const& dS_L_dp_cap,
                      SolidCompressibilityData const& beta_p_SR,
                      FC2aDerivativeData& dfC_2a) const
{
    double const S_L = S_L_data.S_L;
    double const S_G = 1. - S_L;

    double const drho_C_FR_dp_GR =
        /*(dS_G_dp_GR = 0) * constituent_density_data.rho_C_GR +*/
        S_G * phase_transition_data.drho_C_GR_dp_GR +
        /*(dS_L_dp_GR = 0) * constituent_density_data.rho_C_LR +*/
        S_L * phase_transition_data.drho_C_LR_dp_GR;

    dfC_2a.dp_GR = -porosity_data.phi * phase_transition_data.drho_C_GR_dp_GR -
                   drho_C_FR_dp_GR * pCap() *
                       (biot_data() - porosity_data.phi) * beta_p_SR();

    double const dS_G_dp_cap = -dS_L_dp_cap();
    double const rho_C_FR = S_G * constituent_density_data.rho_C_GR +
                            S_L * constituent_density_data.rho_C_LR;

    // TODO (naumov) Extend for partially saturated media.
    constexpr double drho_C_GR_dp_cap = 0;

    double const drho_C_FR_dp_cap =
        dS_G_dp_cap * constituent_density_data.rho_C_GR +
        S_G * drho_C_GR_dp_cap +
        dS_L_dp_cap() * constituent_density_data.rho_C_LR -
        S_L * phase_transition_data.drho_C_LR_dp_LR;

    dfC_2a.dp_cap =
        porosity_data.phi *
            (-phase_transition_data.drho_C_LR_dp_LR - drho_C_GR_dp_cap) -
        drho_C_FR_dp_cap * pCap() * (biot_data() - porosity_data.phi) *
            beta_p_SR() +
        rho_C_FR * (biot_data() - porosity_data.phi) * beta_p_SR();

    double const drho_C_FR_dT = S_G * phase_transition_data.drho_C_GR_dT +
                                S_L * phase_transition_data.drho_C_LR_dT;
    dfC_2a.dT = porosity_d_data.dphi_dT * (constituent_density_data.rho_C_LR -
                                           constituent_density_data.rho_C_GR) +
                porosity_data.phi * (phase_transition_data.drho_C_LR_dT -
                                     phase_transition_data.drho_C_GR_dT) -
                drho_C_FR_dT * pCap() * (biot_data() - porosity_data.phi) *
                    beta_p_SR() +
                rho_C_FR * pCap() * porosity_d_data.dphi_dT * beta_p_SR();
}

void FC3aModel::eval(
    double const dt,
    ConstituentDensityData const& constituent_density_data,
    PrevState<ConstituentDensityData> const& constituent_density_data_prev,
    SaturationData const& S_L_data,
    FC3aData& fC_3a) const
{
    if (dt == 0.)
    {
        fC_3a.a = 0;
        return;
    }

    double const rho_C_GR_dot = (constituent_density_data.rho_C_GR -
                                 constituent_density_data_prev->rho_C_GR) /
                                dt;
    double const rho_C_LR_dot = (constituent_density_data.rho_C_LR -
                                 constituent_density_data_prev->rho_C_LR) /
                                dt;
    auto const S_L = S_L_data.S_L;
    auto const S_G = 1. - S_L;
    fC_3a.a = S_G * rho_C_GR_dot + S_L * rho_C_LR_dot;
}

void FC3aModel::dEval(
    double const dt,
    ConstituentDensityData const& constituent_density_data,
    PrevState<ConstituentDensityData> const& constituent_density_data_prev,
    PhaseTransitionData const& phase_transition_data,
    SaturationData const& S_L_data,
    SaturationDataDeriv const& dS_L_dp_cap,
    FC3aDerivativeData& dfC_3a) const
{
    if (dt == 0.)
    {
        dfC_3a.dp_GR = 0.;
        dfC_3a.dp_cap = 0.;
        dfC_3a.dT = 0.;
        return;
    }
    double const rho_C_GR_dot = (constituent_density_data.rho_C_GR -
                                 constituent_density_data_prev->rho_C_GR) /
                                dt;
    double const rho_C_LR_dot = (constituent_density_data.rho_C_LR -
                                 constituent_density_data_prev->rho_C_LR) /
                                dt;

    auto const S_L = S_L_data.S_L;
    auto const S_G = 1. - S_L;
    dfC_3a.dp_GR =
        /*(dS_G_dp_GR = 0) * rho_C_GR_dot +*/
        S_G * phase_transition_data.drho_C_GR_dp_GR / dt +
        /*(dS_L_dp_GR = 0) * rho_C_LR_dot +*/
        S_L * phase_transition_data.drho_C_LR_dp_GR / dt;

    double const dS_G_dp_cap = -dS_L_dp_cap();
    // TODO (naumov) Extend for partially saturated media.
    constexpr double drho_C_GR_dp_cap = 0;

    dfC_3a.dp_cap = dS_G_dp_cap * rho_C_GR_dot + S_G * drho_C_GR_dp_cap / dt +
                    dS_L_dp_cap() * rho_C_LR_dot -
                    S_L * phase_transition_data.drho_C_LR_dp_LR / dt;

    dfC_3a.dT = S_G * phase_transition_data.drho_C_GR_dT / dt +
                S_L * phase_transition_data.drho_C_LR_dT / dt;
}

template <int DisplacementDim>
void FC4LCpGModel<DisplacementDim>::eval(
    AdvectionData<DisplacementDim> const& advection_data,
    FluidDensityData const& fluid_density_data,
    PhaseTransitionData const& phase_transition_data,
    PorosityData const& porosity_data,
    SaturationData const& S_L_data,
    FC4LCpGData<DisplacementDim>& fC_4_LCpG) const
{
    GlobalDimMatrix<DisplacementDim> const advection_C =
        advection_data.advection_C_G + advection_data.advection_C_L;

    double const sD_G = phase_transition_data.diffusion_coefficient_vapour;
    double const sD_L = phase_transition_data.diffusion_coefficient_solute;

    double const phi_G = (1 - S_L_data.S_L) * porosity_data.phi;
    double const phi_L = S_L_data.S_L * porosity_data.phi;

    double const diffusion_CGpGR = -phi_G * fluid_density_data.rho_GR * sD_G *
                                   phase_transition_data.dxmWG_dpGR;
    double const diffusion_CLpGR = -phi_L * fluid_density_data.rho_LR * sD_L *
                                   phase_transition_data.dxmWL_dpGR;

    double const diffusion_C_pGR = diffusion_CGpGR + diffusion_CLpGR;

    auto const I =
        Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Identity();
    fC_4_LCpG.L.noalias() = diffusion_C_pGR * I + advection_C;
}

template <int DisplacementDim>
void FC4LCpGModel<DisplacementDim>::dEval(
    PermeabilityData<DisplacementDim> const& permeability_data,
    ViscosityData const& viscosity_data,
    PhaseTransitionData const& phase_transition_data,
    AdvectionDerivativeData<DisplacementDim> const& advection_d_data,
    FC4LCpGDerivativeData<DisplacementDim>& dfC_4_LCpG) const
{
    dfC_4_LCpG.dp_GR = advection_d_data.dadvection_C_dp_GR
        // + ddiffusion_C_p_dp_GR TODO (naumov)
        ;

    dfC_4_LCpG.dp_cap = advection_d_data.dadvection_C_dp_cap
        // + ddiffusion_C_p_dp_cap TODO (naumov)
        ;

    GlobalDimMatrix<DisplacementDim> const k_over_mu_G =
        permeability_data.Ki * permeability_data.k_rel_G / viscosity_data.mu_GR;
    GlobalDimMatrix<DisplacementDim> const k_over_mu_L =
        permeability_data.Ki * permeability_data.k_rel_L / viscosity_data.mu_LR;

    dfC_4_LCpG.dT = phase_transition_data.drho_C_GR_dT * k_over_mu_G +
                    phase_transition_data.drho_C_LR_dT * k_over_mu_L
        // + ddiffusion_C_p_dT TODO (naumov)
        ;
}

template struct FC4LCpGModel<2>;
template struct FC4LCpGModel<3>;

template <int DisplacementDim>
void FC4LCpCModel<DisplacementDim>::eval(
    AdvectionData<DisplacementDim> const& advection_data,
    FluidDensityData const& fluid_density_data,
    PhaseTransitionData const& phase_transition_data,
    PorosityData const& porosity_data,
    SaturationData const& S_L_data,
    FC4LCpCData<DisplacementDim>& fC_4_LCpC) const
{
    double const sD_G = phase_transition_data.diffusion_coefficient_vapour;
    double const sD_L = phase_transition_data.diffusion_coefficient_solute;

    double const phi_G = (1 - S_L_data.S_L) * porosity_data.phi;
    double const phi_L = S_L_data.S_L * porosity_data.phi;

    double const diffusion_CGpCap = -phi_G * fluid_density_data.rho_GR * sD_G *
                                    phase_transition_data.dxmWG_dpCap;
    double const diffusion_CLpCap = -phi_L * fluid_density_data.rho_LR * sD_L *
                                    phase_transition_data.dxmWL_dpCap;

    double const diffusion_C_pCap = diffusion_CGpCap + diffusion_CLpCap;

    auto const I =
        Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Identity();

    fC_4_LCpC.L.noalias() = diffusion_C_pCap * I - advection_data.advection_C_L;
}

template <int DisplacementDim>
void FC4LCpCModel<DisplacementDim>::dEval(
    ConstituentDensityData const& constituent_density_data,
    PermeabilityData<DisplacementDim> const& permeability_data,
    PhaseTransitionData const& phase_transition_data,
    SaturationDataDeriv const& dS_L_dp_cap,
    ViscosityData const& viscosity_data,
    FC4LCpCDerivativeData<DisplacementDim>& dfC_4_LCpC) const
{
    ////// Diffusion Part /////
    // TODO (naumov) d(diffusion_C_pCap)/dX for dxmW*/d* != 0

    ////// Advection part /////
    GlobalDimMatrix<DisplacementDim> const k_over_mu_L =
        permeability_data.Ki * permeability_data.k_rel_L / viscosity_data.mu_LR;

    dfC_4_LCpC.dp_GR = phase_transition_data.drho_C_LR_dp_GR * k_over_mu_L
        //+ rhoCLR * (dk_over_mu_L_dp_GR = 0)
        ;

    auto const dk_over_mu_L_dp_cap = permeability_data.Ki *
                                     permeability_data.dk_rel_L_dS_L *
                                     dS_L_dp_cap() / viscosity_data.mu_LR;

    dfC_4_LCpC.dp_cap = -phase_transition_data.drho_C_LR_dp_LR * k_over_mu_L +
                        constituent_density_data.rho_C_LR * dk_over_mu_L_dp_cap;

    dfC_4_LCpC.dT = phase_transition_data.drho_W_LR_dT * k_over_mu_L
        //+ rhoWLR * (dk_over_mu_L_dT != 0 TODO for mu_L(T))
        ;
}

template struct FC4LCpCModel<2>;
template struct FC4LCpCModel<3>;

template <int DisplacementDim>
void FC4LCTModel<DisplacementDim>::eval(
    FluidDensityData const& fluid_density_data,
    PhaseTransitionData const& phase_transition_data,
    PorosityData const& porosity_data,
    SaturationData const& S_L_data,
    FC4LCTData<DisplacementDim>& fC_4_LCT) const
{
    double const sD_G = phase_transition_data.diffusion_coefficient_vapour;
    double const sD_L = phase_transition_data.diffusion_coefficient_solute;

    double const phi_G = (1 - S_L_data.S_L) * porosity_data.phi;
    double const phi_L = S_L_data.S_L * porosity_data.phi;

    double const diffusion_C_G_T = -phi_G * fluid_density_data.rho_GR * sD_G *
                                   phase_transition_data.dxmWG_dT;
    double const diffusion_C_L_T = -phi_L * fluid_density_data.rho_LR * sD_L *
                                   phase_transition_data.dxmWL_dT;

    double const diffusion_C_T = diffusion_C_G_T + diffusion_C_L_T;

    auto const I =
        Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Identity();

    fC_4_LCT.L.noalias() = diffusion_C_T * I;
}

template struct FC4LCTModel<2>;
template struct FC4LCTModel<3>;

void FC4MCpGModel::eval(BiotData const& biot_data,
                        ConstituentDensityData const& constituent_density_data,
                        PorosityData const& porosity_data,
                        SaturationData const& S_L_data,
                        SolidCompressibilityData const& beta_p_SR,
                        FC4MCpGData& fC_4_MCpG) const
{
    auto const S_L = S_L_data.S_L;
    auto const S_G = 1. - S_L;
    double const rho_C_FR = S_G * constituent_density_data.rho_C_GR +
                            S_L * constituent_density_data.rho_C_LR;

    fC_4_MCpG.m = rho_C_FR * (biot_data() - porosity_data.phi) * beta_p_SR();
}

void FC4MCpGModel::dEval(BiotData const& biot_data,
                         ConstituentDensityData const& constituent_density_data,
                         PhaseTransitionData const& phase_transition_data,
                         PorosityData const& porosity_data,
                         PorosityDerivativeData const& porosity_d_data,
                         SaturationData const& S_L_data,
                         SolidCompressibilityData const& beta_p_SR,
                         FC4MCpGDerivativeData& dfC_4_MCpG) const
{
    auto const S_L = S_L_data.S_L;
    auto const S_G = 1. - S_L;
    double const rho_C_FR = S_G * constituent_density_data.rho_C_GR +
                            S_L * constituent_density_data.rho_C_LR;

    double const drho_C_FR_dp_GR =
        /*(dS_G_dp_GR = 0) * constituent_density_data.rho_C_GR +*/
        S_G * phase_transition_data.drho_C_GR_dp_GR +
        /*(dS_L_dp_GR = 0) * constituent_density_data.rho_C_LR +*/
        S_L * phase_transition_data.drho_C_LR_dp_GR;

    dfC_4_MCpG.dp_GR =
        drho_C_FR_dp_GR * (biot_data() - porosity_data.phi) * beta_p_SR();

    double const drho_C_FR_dT = S_G * phase_transition_data.drho_C_GR_dT +
                                S_L * phase_transition_data.drho_C_LR_dT;
    dfC_4_MCpG.dT =
        drho_C_FR_dT * (biot_data() - porosity_data.phi) * beta_p_SR() -
        rho_C_FR * porosity_d_data.dphi_dT * beta_p_SR();
}

void FC4MCpCModel::eval(BiotData const& biot_data,
                        CapillaryPressureData const pCap,
                        ConstituentDensityData const& constituent_density_data,
                        PorosityData const& porosity_data,
                        PrevState<SaturationData> const& S_L_data_prev,
                        SaturationData const& S_L_data,
                        SolidCompressibilityData const& beta_p_SR,
                        FC4MCpCData& fC_4_MCpC) const
{
    auto const S_L = S_L_data.S_L;
    auto const S_G = 1. - S_L;
    double const rho_C_FR = S_G * constituent_density_data.rho_C_GR +
                            S_L * constituent_density_data.rho_C_LR;

    fC_4_MCpC.m =
        -rho_C_FR * (biot_data() - porosity_data.phi) * beta_p_SR() * S_L;

    fC_4_MCpC.ml =
        (porosity_data.phi * (constituent_density_data.rho_C_LR -
                              constituent_density_data.rho_C_GR) -
         rho_C_FR * pCap() * (biot_data() - porosity_data.phi) * beta_p_SR()) *
        (S_L - S_L_data_prev->S_L);
}

template <int DisplacementDim>
void FC4MCTModel<DisplacementDim>::eval(
    BiotData const& biot_data,
    ConstituentDensityData const& constituent_density_data,
    PorosityData const& porosity_data,
    SaturationData const& S_L_data,
    SolidThermalExpansionData<DisplacementDim> const& s_therm_exp_data,
    FC4MCTData& fC_4_MCT) const
{
    auto const S_L = S_L_data.S_L;
    auto const S_G = 1. - S_L;
    double const rho_C_FR = S_G * constituent_density_data.rho_C_GR +
                            S_L * constituent_density_data.rho_C_LR;

    fC_4_MCT.m = -rho_C_FR * (biot_data() - porosity_data.phi) *
                 s_therm_exp_data.beta_T_SR;
}

template <int DisplacementDim>
void FC4MCTModel<DisplacementDim>::dEval(
    BiotData const& biot_data,
    [[maybe_unused]] ConstituentDensityData const& constituent_density_data,
    PhaseTransitionData const& phase_transition_data,
    PorosityData const& porosity_data,
    [[maybe_unused]] PorosityDerivativeData const& porosity_d_data,
    SaturationData const& S_L_data,
    SolidThermalExpansionData<DisplacementDim> const& s_therm_exp_data,
    FC4MCTDerivativeData& dfC_4_MCT) const
{
    auto const S_L = S_L_data.S_L;
    auto const S_G = 1. - S_L;
#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
    double const rho_C_FR = S_G * constituent_density_data.rho_C_GR +
                            S_L * constituent_density_data.rho_C_LR;
#endif

    double const drho_C_FR_dT = S_G * phase_transition_data.drho_C_GR_dT +
                                S_L * phase_transition_data.drho_C_LR_dT;

    dfC_4_MCT.dT = drho_C_FR_dT * (biot_data() - porosity_data.phi) *
                       s_therm_exp_data.beta_T_SR
#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
                   + rho_C_FR * (biot_data() - porosity_d_data.dphi_dT) *
                         s_therm_exp_data.beta_T_SR
#endif
        ;
}

template struct FC4MCTModel<2>;
template struct FC4MCTModel<3>;

void FC4MCuModel::eval(BiotData const& biot_data,
                       ConstituentDensityData const& constituent_density_data,
                       SaturationData const& S_L_data,
                       FC4MCuData& fC_4_MCu) const
{
    auto const S_L = S_L_data.S_L;
    auto const S_G = 1. - S_L;
    double const rho_C_FR = S_G * constituent_density_data.rho_C_GR +
                            S_L * constituent_density_data.rho_C_LR;

    fC_4_MCu.m = rho_C_FR * biot_data();
}

void FC4MCuModel::dEval(BiotData const& biot_data,
                        PhaseTransitionData const& phase_transition_data,
                        SaturationData const& S_L_data,
                        FC4MCuDerivativeData& dfC_4_MCu) const
{
    auto const S_L = S_L_data.S_L;
    auto const S_G = 1. - S_L;
    double const drho_C_FR_dT = S_G * phase_transition_data.drho_C_GR_dT +
                                S_L * phase_transition_data.drho_C_LR_dT;
    dfC_4_MCu.dT = drho_C_FR_dT * biot_data();
}
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
