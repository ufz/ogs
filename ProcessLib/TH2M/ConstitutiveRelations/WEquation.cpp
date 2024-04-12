/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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
        rho_W_FR * pCap() * (biot_data() - porosity_data.phi) * beta_p_SR();
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
                  drho_C_FR_dp_GR * pCap() * (biot_data() - porosity_data.phi) *
                      beta_p_SR();

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
        drho_W_FR_dp_cap * pCap() * (biot_data() - porosity_data.phi) *
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
        drho_W_FR_dT * pCap() * (biot_data() - porosity_data.phi) *
            beta_p_SR() -
        rho_W_FR * pCap() * porosity_d_data.dphi_dT * beta_p_SR();

    dfW_2.dT = dfW_2a_dT - dfW_2b_dT;
}
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
