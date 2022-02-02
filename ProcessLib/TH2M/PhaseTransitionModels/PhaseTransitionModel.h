/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <map>

#include "MaterialLib/MPL/Medium.h"

namespace ProcessLib
{
namespace TH2M
{
struct PhaseTransitionModelVariables
{
    // gas phase density
    double rhoGR = 0.;
    double rhoCGR = 0.;
    double rhoWGR = 0.;

    double drho_GR_dp_GR = 0.;
    double drho_GR_dp_cap = 0.;
    double drho_GR_dT = 0.;

    double drho_C_GR_dp_GR = 0.;
    double drho_C_GR_dp_cap = 0.;
    double drho_C_GR_dT = 0.;

    double drho_W_GR_dp_GR = 0.;
    double drho_W_GR_dp_cap = 0.;
    double drho_W_GR_dT = 0.;

    // liquid phase density
    double rhoLR = 0.;
    double rhoWLR = 0.;
    double rhoCLR = 0.;

    double drho_LR_dp_GR = 0.;
    double drho_LR_dp_cap = 0.;
    double drho_LR_dT = 0.;
    double drho_LR_dp_LR = 0.;

    double drho_C_LR_dp_GR = 0.;
    double drho_C_LR_dp_cap = 0.;
    double drho_C_LR_dT = 0.;
    double drho_C_LR_dp_LR = 0.;

    double drho_W_LR_dp_GR = 0.;
    double drho_W_LR_dp_cap = 0.;
    double drho_W_LR_dT = 0.;
    double drho_W_LR_dp_LR = 0.;

    // water partial pressure in gas phase
    double pWGR = 0;

    // constituent mass and molar fractions
    double xnWG = 0.;
    double xmWG = 0.;
    double xnWL = 0.;
    double xmWL = 0.;

    // mass fraction derivatives
    double dxmWG_dpGR = 0.;
    double dxmWG_dpCap = 0.;
    double dxmWG_dT = 0.;

    double dxmWL_dpGR = 0.;
    double dxmWL_dpCap = 0.;
    double dxmWL_dpLR = 0.;
    double dxmWL_dT = 0.;

    // viscosities
    double muGR = 0.;
    double muLR = 0.;

    double diffusion_coefficient_vapour = 0.;
    double diffusion_coefficient_solvate = 0.;

    // specific enthalpies
    double hG = 0;
    double hCG = 0;
    double hWG = 0;
    double hL = 0;

    double dh_G_dT = 0;
    double dh_L_dT = 0;

    // specific inner energies
    double uG = 0;
    double uL = 0;

    double du_G_dT = 0;
    double du_L_dT = 0;
    double du_G_dp_GR = 0;
    double du_L_dp_GR = 0;
    double du_L_dp_cap = 0;
};

struct PhaseTransitionModel
{
    explicit PhaseTransitionModel(
        std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const&
            media)
    {
        DBUG("Create phase transition models...");

        // check for minimum requirement definitions in media object
        std::array const required_gas_properties = {
            MaterialPropertyLib::viscosity, MaterialPropertyLib::density};
        std::array const required_liquid_properties = {
            MaterialPropertyLib::viscosity, MaterialPropertyLib::density};

        for (auto const& m : media)
        {
            checkRequiredProperties(m.second->phase("Gas"),
                                    required_gas_properties);
            checkRequiredProperties(m.second->phase("AqueousLiquid"),
                                    required_liquid_properties);
        }
    }

    virtual ~PhaseTransitionModel() = default;

    virtual void computeConstitutiveVariables(
        const MaterialPropertyLib::Medium* medium,
        MaterialPropertyLib::VariableArray variables,
        ParameterLib::SpatialPosition pos, double const t, double const dt)
    {
        cv = updateConstitutiveVariables(cv, medium, variables, pos, t, dt);
    }

    virtual PhaseTransitionModelVariables updateConstitutiveVariables(
        PhaseTransitionModelVariables const& phase_transition_model_variables,
        const MaterialPropertyLib::Medium* medium,
        MaterialPropertyLib::VariableArray variables,
        ParameterLib::SpatialPosition pos, double const t,
        double const dt) const = 0;

    // constitutive variables
    PhaseTransitionModelVariables cv;
};

int numberOfComponents(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media,
    std::string phase_name);

int findComponentIndex(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media,
    std::string phase_name, MaterialPropertyLib::PropertyType property_type);

}  // namespace TH2M
}  // namespace ProcessLib
