/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "WEquation.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
void FW1Model<DisplacementDim>::eval(
    AdvectionData<DisplacementDim> const& advection_data,
    FluidDensityData const& fluid_density_data,
    FW1Data<DisplacementDim>& fW_1) const
{
    fW_1.A = advection_data.advection_W_G * fluid_density_data.rho_GR +
             advection_data.advection_W_L * fluid_density_data.rho_LR;
}

template struct FW1Model<2>;
template struct FW1Model<3>;

void FW2Model::eval(BiotData const biot_data,
                    CapillaryPressureData const pCap,
                    ConstituentDensityData const& constituent_density_data,
                    PorosityData const& porosity_data,
                    PureLiquidDensityData const& rho_W_LR,
                    SaturationData const& S_L_data,
                    SolidCompressibilityData const beta_p_SR,
                    FW2Data& fW_2) const
{
    auto const S_L = S_L_data.S_L;
    auto const S_G = 1. - S_L;
    double const rho_W_FR =
        S_G * constituent_density_data.rho_W_GR + S_L * rho_W_LR();

    fW_2.a =
        porosity_data.phi * (rho_W_LR() - constituent_density_data.rho_W_GR) -
        rho_W_FR * pCap.pCap * (biot_data() - porosity_data.phi) * beta_p_SR();
}

void FW2Model::dEval(BiotData const& biot_data,
                     CapillaryPressureData const pCap,
                     ConstituentDensityData const& constituent_density_data,
                     PhaseTransitionData const& phase_transition_data,
                     PorosityData const& porosity_data,
                     PorosityDerivativeData const& porosity_d_data,
                     PureLiquidDensityData const& rho_W_LR,
                     SaturationData const& S_L_data,
                     SaturationDataDeriv const& dS_L_dp_cap,
                     SolidCompressibilityData const& beta_p_SR,
                     FW2DerivativeData& dfW_2) const
{
    double const S_L = S_L_data.S_L;
    double const S_G = 1. - S_L;

    double const drho_C_FR_dp_GR =
        /*(dS_G_dp_GR = 0) * constituent_density_data.rho_C_GR +*/
        S_G * phase_transition_data.drho_C_GR_dp_GR +
        /*(dS_L_dp_GR = 0) * constituent_density_data.rho_C_LR +*/
        S_L * phase_transition_data.drho_C_LR_dp_GR;

    dfW_2.dp_GR = -porosity_data.phi * phase_transition_data.drho_C_GR_dp_GR -
                  drho_C_FR_dp_GR * pCap.pCap *
                      (biot_data() - porosity_data.phi) * beta_p_SR();

    double const dfW_2a_dp_cap =
        porosity_data.phi * (-phase_transition_data.drho_W_LR_dp_LR -
                             phase_transition_data.drho_W_GR_dp_cap);
    double const rho_W_FR =
        S_G * constituent_density_data.rho_W_GR + S_L * rho_W_LR();

    double const drho_W_FR_dp_cap =
        -dS_L_dp_cap() * constituent_density_data.rho_W_GR +
        S_G * phase_transition_data.drho_W_GR_dp_cap +
        dS_L_dp_cap() * rho_W_LR() -
        S_L * phase_transition_data.drho_W_LR_dp_LR;

    double const dfW_2b_dp_cap =
        drho_W_FR_dp_cap * pCap.pCap * (biot_data() - porosity_data.phi) *
            beta_p_SR() +
        rho_W_FR * (biot_data() - porosity_data.phi) * beta_p_SR();

    dfW_2.dp_cap = dfW_2a_dp_cap - dfW_2b_dp_cap;

    double const drho_W_FR_dT = S_G * phase_transition_data.drho_W_GR_dT +
                                S_L * phase_transition_data.drho_W_LR_dT;

    double const dfW_2a_dT =
        porosity_d_data.dphi_dT *
            (rho_W_LR() - constituent_density_data.rho_W_GR) +
        porosity_data.phi * (phase_transition_data.drho_W_LR_dT -
                             phase_transition_data.drho_W_GR_dT);
    double const dfW_2b_dT =
        drho_W_FR_dT * pCap.pCap * (biot_data() - porosity_data.phi) *
            beta_p_SR() -
        rho_W_FR * pCap.pCap * porosity_d_data.dphi_dT * beta_p_SR();

    dfW_2.dT = dfW_2a_dT - dfW_2b_dT;
}

void FW3aModel::eval(
    double const dt,
    ConstituentDensityData const& constituent_density_data,
    PrevState<ConstituentDensityData> const& constituent_density_data_prev,
    PrevState<PureLiquidDensityData> const& rho_W_LR_prev,
    PureLiquidDensityData const& rho_W_LR,
    SaturationData const& S_L_data,
    FW3aData& fW_3a) const
{
    if (dt == 0.)
    {
        fW_3a.a = 0;
        return;
    }

    double const rho_W_GR_dot = (constituent_density_data.rho_W_GR -
                                 constituent_density_data_prev->rho_W_GR) /
                                dt;
    double const rho_W_LR_dot = (rho_W_LR() - **rho_W_LR_prev) / dt;
    auto const S_L = S_L_data.S_L;
    auto const S_G = 1. - S_L;
    fW_3a.a = S_G * rho_W_GR_dot + S_L * rho_W_LR_dot;
}

void FW3aModel::dEval(
    double const dt,
    ConstituentDensityData const& constituent_density_data,
    PhaseTransitionData const& phase_transition_data,
    PrevState<ConstituentDensityData> const& constituent_density_data_prev,
    PrevState<PureLiquidDensityData> const& rho_W_LR_prev,
    PureLiquidDensityData const& rho_W_LR,
    SaturationData const& S_L_data,
    SaturationDataDeriv const& dS_L_dp_cap,
    FW3aDerivativeData& dfW_3a) const
{
    if (dt == 0.)
    {
        dfW_3a.dp_GR = 0.;
        dfW_3a.dp_cap = 0.;
        dfW_3a.dT = 0.;
        return;
    }

    auto const S_L = S_L_data.S_L;
    auto const S_G = 1. - S_L;

    double const rho_W_GR_dot = (constituent_density_data.rho_W_GR -
                                 constituent_density_data_prev->rho_W_GR) /
                                dt;
    double const rho_W_LR_dot = (rho_W_LR() - **rho_W_LR_prev) / dt;

    dfW_3a.dp_GR = /*(ds_G_dp_GR = 0) * rho_W_GR_dot +*/
        S_G * phase_transition_data.drho_W_GR_dp_GR / dt +
        /*(ds_L_dp_GR = 0) * rho_W_LR_dot +*/
        S_L * phase_transition_data.drho_W_LR_dp_GR / dt;

    dfW_3a.dp_cap = -dS_L_dp_cap() * rho_W_GR_dot +
                    S_G * phase_transition_data.drho_W_GR_dp_cap / dt +
                    dS_L_dp_cap() * rho_W_LR_dot -
                    S_L * phase_transition_data.drho_W_LR_dp_LR / dt;
    dfW_3a.dT = S_G * phase_transition_data.drho_W_GR_dT / dt +
                S_L * phase_transition_data.drho_W_LR_dT / dt;
}

template <int DisplacementDim>
void FW4LWpGModel<DisplacementDim>::eval(
    AdvectionData<DisplacementDim> const& advection_data,
    FluidDensityData const& fluid_density_data,
    PhaseTransitionData const& phase_transition_data,
    PorosityData const& porosity_data,
    SaturationData const& S_L_data,
    FW4LWpGData<DisplacementDim>& fW_4_LWpG) const
{
    GlobalDimMatrix<DisplacementDim> const advection_W =
        advection_data.advection_W_G + advection_data.advection_W_L;

    double const sD_G = phase_transition_data.diffusion_coefficient_vapour;
    double const sD_L = phase_transition_data.diffusion_coefficient_solute;

    double const phi_G = (1 - S_L_data.S_L) * porosity_data.phi;
    double const phi_L = S_L_data.S_L * porosity_data.phi;

    double const diffusion_WGpGR = phi_G * fluid_density_data.rho_GR * sD_G *
                                   phase_transition_data.dxmWG_dpGR;
    double const diffusion_WLpGR = phi_L * fluid_density_data.rho_LR * sD_L *
                                   phase_transition_data.dxmWL_dpGR;
    double const diffusion_W_pGR = diffusion_WGpGR + diffusion_WLpGR;

    auto const I =
        Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Identity();
    fW_4_LWpG.L.noalias() = diffusion_W_pGR * I + advection_W;
}

template <int DisplacementDim>
void FW4LWpGModel<DisplacementDim>::dEval(
    ConstituentDensityData const& constituent_density_data,
    PermeabilityData<DisplacementDim> const& permeability_data,
    PhaseTransitionData const& phase_transition_data,
    PureLiquidDensityData const& rho_W_LR,
    SaturationDataDeriv const& dS_L_dp_cap,
    ViscosityData const& viscosity_data,
    FW4LWpGDerivativeData<DisplacementDim>& dfW_4_LWpG) const
{
    ////// Diffusion Part /////
    // TODO (naumov) d(diffusion_W_p)/dX for dxmW*/d* != 0

    auto const k_over_mu_G =
        permeability_data.Ki * permeability_data.k_rel_G / viscosity_data.mu_GR;
    auto const k_over_mu_L =
        permeability_data.Ki * permeability_data.k_rel_L / viscosity_data.mu_LR;

    // dk_over_mu_G_dp_GR = ip_out.permeability_data.Ki *
    //                      ip_out.permeability_data.dk_rel_G_dS_L *
    //                      (ds_L_dp_GR = 0) /
    //                      ip_cv.viscosity_data.mu_GR = 0;
    // dk_over_mu_L_dp_GR = ip_out.permeability_data.Ki *
    //                      ip_out.permeability_data.dk_rel_L_dS_L *
    //                      (ds_L_dp_GR = 0) /
    //                      ip_cv.viscosity_data.mu_LR = 0;
    auto const dk_over_mu_G_dp_cap = permeability_data.Ki *
                                     permeability_data.dk_rel_G_dS_L *
                                     dS_L_dp_cap() / viscosity_data.mu_GR;

    auto const dk_over_mu_L_dp_cap = permeability_data.Ki *
                                     permeability_data.dk_rel_L_dS_L *
                                     dS_L_dp_cap() / viscosity_data.mu_LR;

    dfW_4_LWpG.dp_GR = phase_transition_data.drho_W_GR_dp_GR * k_over_mu_G
                       // + rhoWGR * (dk_over_mu_G_dp_GR = 0)
                       + phase_transition_data.drho_W_LR_dp_GR * k_over_mu_L
        // + rhoWLR * (dk_over_mu_L_dp_GR = 0)
        ;

    dfW_4_LWpG.dp_cap =
        phase_transition_data.drho_W_GR_dp_cap * k_over_mu_G +
        constituent_density_data.rho_W_GR * dk_over_mu_G_dp_cap +
        -phase_transition_data.drho_W_LR_dp_LR * k_over_mu_L +
        rho_W_LR() * dk_over_mu_L_dp_cap;

    dfW_4_LWpG.dT = phase_transition_data.drho_W_GR_dT * k_over_mu_G
                    //+ rhoWGR * (dk_over_mu_G_dT != 0 TODO for mu_G(T))
                    + phase_transition_data.drho_W_LR_dT * k_over_mu_L
        //+ rhoWLR * (dk_over_mu_L_dT != 0 TODO for mu_G(T))
        ;
}

template struct FW4LWpGModel<2>;
template struct FW4LWpGModel<3>;

template <int DisplacementDim>
void FW4LWpCModel<DisplacementDim>::eval(
    AdvectionData<DisplacementDim> const& advection_data,
    FluidDensityData const& fluid_density_data,
    PhaseTransitionData const& phase_transition_data,
    PorosityData const& porosity_data,
    SaturationData const& S_L_data,
    FW4LWpCData<DisplacementDim>& fW_4_LWpC) const
{
    double const sD_G = phase_transition_data.diffusion_coefficient_vapour;
    double const sD_L = phase_transition_data.diffusion_coefficient_solute;

    double const phi_G = (1 - S_L_data.S_L) * porosity_data.phi;
    double const phi_L = S_L_data.S_L * porosity_data.phi;

    double const diffusion_WGpCap = phi_G * fluid_density_data.rho_GR * sD_G *
                                    phase_transition_data.dxmWG_dpCap;
    double const diffusion_WLpCap = phi_L * fluid_density_data.rho_LR * sD_L *
                                    phase_transition_data.dxmWL_dpCap;

    double const diffusion_W_pCap = diffusion_WGpCap + diffusion_WLpCap;

    auto const I =
        Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Identity();

    fW_4_LWpC.L.noalias() = diffusion_W_pCap * I - advection_data.advection_W_L;
}

template <int DisplacementDim>
void FW4LWpCModel<DisplacementDim>::dEval(
    AdvectionData<DisplacementDim> const& advection_data,
    FluidDensityData const& fluid_density_data,
    PermeabilityData<DisplacementDim> const& permeability_data,
    PhaseTransitionData const& phase_transition_data,
    PorosityData const& porosity_data,
    PureLiquidDensityData const& rho_W_LR,
    SaturationData const& S_L_data,
    SaturationDataDeriv const& dS_L_dp_cap,
    ViscosityData const& viscosity_data,
    FW4LWpCDerivativeData<DisplacementDim>& dfW_4_LWpC) const
{
    ////// Diffusion Part /////
    // TODO (naumov) d(diffusion_W_pCap)/dX for dxmW*/d* != 0

    ////// Advection part /////
    GlobalDimMatrix<DisplacementDim> const k_over_mu_L =
        permeability_data.Ki * permeability_data.k_rel_L / viscosity_data.mu_LR;

    dfW_4_LWpC.dp_GR = phase_transition_data.drho_W_LR_dp_GR * k_over_mu_L
        //+ rhoWLR * (dk_over_mu_L_dp_GR = 0)
        ;

    double const sD_G = phase_transition_data.diffusion_coefficient_vapour;
    double const sD_L = phase_transition_data.diffusion_coefficient_solute;

    double const phi_G = (1 - S_L_data.S_L) * porosity_data.phi;
    double const phi_L = S_L_data.S_L * porosity_data.phi;

    double const diffusion_WGpCap = phi_G * fluid_density_data.rho_GR * sD_G *
                                    phase_transition_data.dxmWG_dpCap;
    double const diffusion_WLpCap = phi_L * fluid_density_data.rho_LR * sD_L *
                                    phase_transition_data.dxmWL_dpCap;

    double const diffusion_W_pCap = diffusion_WGpCap + diffusion_WLpCap;

    auto const I =
        Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Identity();

    dfW_4_LWpC.dp_cap = diffusion_W_pCap * I - advection_data.advection_W_L;

    auto const dk_over_mu_L_dp_cap = permeability_data.Ki *
                                     permeability_data.dk_rel_L_dS_L *
                                     dS_L_dp_cap() / viscosity_data.mu_LR;
    dfW_4_LWpC.dp_cap = -phase_transition_data.drho_W_LR_dp_LR * k_over_mu_L +
                        rho_W_LR() * dk_over_mu_L_dp_cap;

    dfW_4_LWpC.dT = phase_transition_data.drho_W_LR_dT * k_over_mu_L
        //+ rhoWLR * (dk_over_mu_L_dT != 0 TODO for mu_L(T))
        ;
}

template struct FW4LWpCModel<2>;
template struct FW4LWpCModel<3>;

template <int DisplacementDim>
void FW4LWTModel<DisplacementDim>::eval(
    FluidDensityData const& fluid_density_data,
    PhaseTransitionData const& phase_transition_data,
    PorosityData const& porosity_data,
    SaturationData const& S_L_data,
    FW4LWTData<DisplacementDim>& fW_4_LWT) const
{
    double const sD_G = phase_transition_data.diffusion_coefficient_vapour;
    double const sD_L = phase_transition_data.diffusion_coefficient_solute;

    double const phi_G = (1 - S_L_data.S_L) * porosity_data.phi;
    double const phi_L = S_L_data.S_L * porosity_data.phi;

    double const diffusion_W_G_T = phi_G * fluid_density_data.rho_GR * sD_G *
                                   phase_transition_data.dxmWG_dT;
    double const diffusion_W_L_T = phi_L * fluid_density_data.rho_LR * sD_L *
                                   phase_transition_data.dxmWL_dT;

    double const diffusion_W_T = diffusion_W_G_T + diffusion_W_L_T;

    auto const I =
        Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Identity();

    fW_4_LWT.L.noalias() = diffusion_W_T * I;
}

template struct FW4LWTModel<2>;
template struct FW4LWTModel<3>;

void FW4MWpGModel::eval(BiotData const& biot_data,
                        ConstituentDensityData const& constituent_density_data,
                        PorosityData const& porosity_data,
                        PureLiquidDensityData const& rho_W_LR,
                        SaturationData const& S_L_data,
                        SolidCompressibilityData const& beta_p_SR,
                        FW4MWpGData& fW_4_MWpG) const
{
    double const S_L = S_L_data.S_L;
    double const S_G = 1. - S_L;

    double const rho_W_FR =
        S_G * constituent_density_data.rho_W_GR + S_L * rho_W_LR();

    fW_4_MWpG.m = rho_W_FR * (biot_data() - porosity_data.phi) * beta_p_SR();
}

void FW4MWpCModel::eval(BiotData const& biot_data,
                        CapillaryPressureData const pCap,
                        ConstituentDensityData const& constituent_density_data,
                        PorosityData const& porosity_data,
                        PrevState<SaturationData> const& S_L_data_prev,
                        PureLiquidDensityData const& rho_W_LR,
                        SaturationData const& S_L_data,
                        SolidCompressibilityData const& beta_p_SR,
                        FW4MWpCData& fW_4_MWpC) const
{
    auto const S_L = S_L_data.S_L;
    auto const S_G = 1. - S_L;
    double const rho_W_FR =
        S_G * constituent_density_data.rho_W_GR + S_L * rho_W_LR();

    fW_4_MWpC.m =
        -rho_W_FR * (biot_data() - porosity_data.phi) * beta_p_SR() * S_L;

    fW_4_MWpC.ml =
        (porosity_data.phi * (rho_W_LR() - constituent_density_data.rho_W_GR) -
         rho_W_FR * pCap.pCap * (biot_data() - porosity_data.phi) *
             beta_p_SR()) *
        (S_L - S_L_data_prev->S_L);
}

template <int DisplacementDim>
void FW4MWTModel<DisplacementDim>::eval(
    BiotData const& biot_data,
    ConstituentDensityData const& constituent_density_data,
    PorosityData const& porosity_data,
    PureLiquidDensityData const& rho_W_LR,
    SaturationData const& S_L_data,
    SolidThermalExpansionData<DisplacementDim> const& s_therm_exp_data,
    FW4MWTData& fW_4_MWT) const
{
    auto const S_L = S_L_data.S_L;
    auto const S_G = 1. - S_L;
    double const rho_W_FR =
        S_G * constituent_density_data.rho_W_GR + S_L * rho_W_LR();

    fW_4_MWT.m = -rho_W_FR * (biot_data() - porosity_data.phi) *
                 s_therm_exp_data.beta_T_SR;
}

template struct FW4MWTModel<2>;
template struct FW4MWTModel<3>;

void FW4MWuModel::eval(BiotData const& biot_data,
                       ConstituentDensityData const& constituent_density_data,
                       PureLiquidDensityData const& rho_W_LR,
                       SaturationData const& S_L_data,
                       FW4MWuData& fW_4_MWu) const
{
    auto const S_L = S_L_data.S_L;
    auto const S_G = 1. - S_L;
    double const rho_W_FR =
        S_G * constituent_density_data.rho_W_GR + S_L * rho_W_LR();

    fW_4_MWu.m = rho_W_FR * biot_data();
}

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
