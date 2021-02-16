/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <map>

#include "MaterialLib/MPL/Medium.h"

namespace ProcessLib
{
namespace TH2M
{
struct PhaseTransitionModels
{
    explicit PhaseTransitionModels(
        std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const&
            media)
    {
        DBUG("Create phase transition models...");

        // check for minimum requirement definitions in media object
        std::array const required_gas_properties = {
            MaterialPropertyLib::viscosity, MaterialPropertyLib::density,
            MaterialPropertyLib::thermal_conductivity};
        std::array const required_liquid_properties = {
            MaterialPropertyLib::viscosity, MaterialPropertyLib::density,
            MaterialPropertyLib::thermal_conductivity};

        for (auto const& m : media)
        {
            checkRequiredProperties(m.second->phase("Gas"),
                                    required_gas_properties);
            checkRequiredProperties(m.second->phase("AqueousLiquid"),
                                    required_liquid_properties);
        }
    }

    virtual ~PhaseTransitionModels() = default;

    virtual void getConstitutiveVariables(
        const MaterialPropertyLib::Medium* medium,
        MaterialPropertyLib::VariableArray variables,
        ParameterLib::SpatialPosition pos, double const t, double const dt) = 0;

    // constitutive variables as public members

    // gas phase density
    double rhoGR = 0.;
    double rhoCGR = 0.;
    double rhoWGR = 0.;

    // liquid phase density
    double rhoLR = 0.;
    double rhoWLR = 0.;
    double rhoCLR = 0.;

    // water partial pressure in gas phase
    double pWGR = 0;

    // constituent mass and molar fractions
    double xnCG = 0.;
    double xnWG = 0.;
    double xmCG = 0.;
    double xmWG = 0.;
    double xmCL = 0.;
    double xmWL = 0.;

    // molar fraction derivatives
    double dxnCG_dpGR = 0.;
    double dxnCG_dpCap = 0.;
    double dxnCG_dT = 0.;

    // mass fraction derivatives
    double dxmCG_dpGR = 0.;
    double dxmWG_dpGR = 0.;
    double dxmCL_dpLR = 0.;
    double dxmWL_dpLR = 0.;
    double dxmCG_dT = 0.;
    double dxmWG_dT = 0.;
    double dxmCL_dT = 0.;
    double dxmWL_dT = 0.;

    // viscosities
    double muGR = 0.;
    double muLR = 0.;
    // thermal conductivities
    double lambdaGR = 0.;
    double lambdaLR = 0.;

    double diffusion_coefficient_vapour = 0.;
    double diffusion_coefficient_solvate = 0.;

    // specific enthalpies
    double hG = 0;
    double hCG = 0;
    double hWG = 0;
    double hL = 0;

    // specific inner energies
    double uG = 0;
    double uL = 0;
};

struct PhaseTransitionNone : PhaseTransitionModels
{
    explicit PhaseTransitionNone(
        std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const&
            media);

    void getConstitutiveVariables(const MaterialPropertyLib::Medium* medium,
                                  MaterialPropertyLib::VariableArray variables,
                                  ParameterLib::SpatialPosition pos,
                                  double const t, double const dt) override;
};

struct PhaseTransitionEvaporation : PhaseTransitionModels
{
    explicit PhaseTransitionEvaporation(
        std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const&
            media);

    void getConstitutiveVariables(const MaterialPropertyLib::Medium* medium,
                                  MaterialPropertyLib::VariableArray variables,
                                  ParameterLib::SpatialPosition pos,
                                  double const t, double const dt) override;

private:
    std::size_t n_components_gas_;
    int gas_phase_vapour_component_index_ = -1;
    int gas_phase_dry_air_component_index_ = -1;
};

// Dissolution only: Amounts of the gas phase can be dissolved into the liquid
// phase. This is realized by defining two components in the liqid phase but
// only one (or zero) components in the gas phase.
struct PhaseTransitionDissolution : PhaseTransitionModels
{
    explicit PhaseTransitionDissolution(
        std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const&
            media);

    void getConstitutiveVariables(const MaterialPropertyLib::Medium* medium,
                                  MaterialPropertyLib::VariableArray variables,
                                  ParameterLib::SpatialPosition pos,
                                  double const t, double const dt) override;
};

// Full phase transition: Gas can dissolve into the liquid phase according to an
// equilibrium and water evaporates into the gas phase according to Dalton's
// law. This is realized by defining two components in each gas and liquid
// phase.
struct PhaseTransitionFull : PhaseTransitionModels
{
    explicit PhaseTransitionFull(
        std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const&
            media);

    void getConstitutiveVariables(const MaterialPropertyLib::Medium* medium,
                                  MaterialPropertyLib::VariableArray variables,
                                  ParameterLib::SpatialPosition pos,
                                  double const t, double const dt) override;
};

}  // namespace TH2M
}  // namespace ProcessLib
