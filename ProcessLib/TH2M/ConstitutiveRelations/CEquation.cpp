/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CEquation.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
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
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
