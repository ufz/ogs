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

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
struct PhaseTransitionData
{
    double drho_GR_dp_GR = nan;
    double drho_GR_dp_cap = nan;
    double drho_GR_dT = nan;

    double drho_C_GR_dp_GR = nan;
    double drho_C_GR_dp_cap = nan;
    double drho_C_GR_dT = nan;

    double drho_W_GR_dp_GR = nan;
    double drho_W_GR_dp_cap = nan;
    double drho_W_GR_dT = nan;

    double drho_LR_dT = nan;
    double drho_LR_dp_LR = nan;

    // TODO (naumov) These three are zero in both models but used in the
    // assembly. Remove them and simplify assembly or correct the expressions in
    // the phase transition models?
    static constexpr double drho_C_LR_dp_GR = 0.;
    static constexpr double drho_C_LR_dT = 0.;
    static constexpr double drho_C_LR_dp_LR = 0.;

    double drho_W_LR_dp_GR = nan;
    double drho_W_LR_dT = nan;
    double drho_W_LR_dp_LR = nan;

    // mass fraction derivatives
    double dxmWG_dpGR = nan;
    double dxmWG_dpCap = nan;
    double dxmWG_dT = nan;

    double dxmWL_dpGR = nan;
    double dxmWL_dpCap = nan;
    double dxmWL_dT = nan;

    double diffusion_coefficient_vapour = nan;
    double diffusion_coefficient_solute = nan;

    // specific enthalpies
    double hCG = nan;
    double hWG = nan;

    double dh_G_dT = nan;
    double dh_L_dT = nan;

    // specific inner energies
    double uG = nan;
    double uL = nan;

    double du_G_dT = nan;
    double du_L_dT = nan;
    double du_G_dp_GR = nan;

    // TODO (naumov) These two are zero in both models but used in the assembly.
    // Remove them and simplify assembly or correct the expressions in the phase
    // transition models?
    static constexpr double du_L_dp_GR = 0;
    static constexpr double du_L_dp_cap = 0;
};

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
