/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "PhaseTransitionNone.h"

#include "MaterialLib/PhysicalConstant.h"

namespace ProcessLib
{
namespace TH2M
{
PhaseTransitionNone::PhaseTransitionNone(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
    : PhaseTransitionModels(media)
{
    DBUG("Create PhaseTransitionNone constitutive model.");

    // check for minimum requirement definitions in media object
    std::array const required_gas_properties = {
        MaterialPropertyLib::specific_heat_capacity,
        MaterialPropertyLib::molar_mass};
    std::array const required_liquid_properties = {
        MaterialPropertyLib::specific_heat_capacity};

    for (auto const& m : media)
    {
        checkRequiredProperties(m.second->phase("Gas"),
                                required_gas_properties);
        checkRequiredProperties(m.second->phase("AqueousLiquid"),
                                required_liquid_properties);
    }
}

void PhaseTransitionNone::computeConstitutiveVariables(
    const MaterialPropertyLib::Medium* medium,
    MaterialPropertyLib::VariableArray variables,
    ParameterLib::SpatialPosition pos, double const t, double const dt)
{
    // primary variables
    auto const pGR = std::get<double>(variables[static_cast<int>(
        MaterialPropertyLib::Variable::phase_pressure)]);
    auto const T = std::get<double>(variables[static_cast<int>(
        MaterialPropertyLib::Variable::temperature)]);

    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& gas_phase = medium->phase("Gas");

    // C-component is only component in the gas phase
    xnCG = 1.;
    xmCG = 1.;

    auto const M =
        gas_phase.property(MaterialPropertyLib::PropertyType::molar_mass)
            .template value<double>(variables, pos, t, dt);

    variables[static_cast<int>(MaterialPropertyLib::Variable::molar_mass)] = M;

    rhoGR = gas_phase.property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(variables, pos, t, dt);
    muGR = gas_phase.property(MaterialPropertyLib::PropertyType::viscosity)
               .template value<double>(variables, pos, t, dt);
    lambdaGR =
        gas_phase
            .property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .template value<double>(variables, pos, t, dt);

    rhoCGR = rhoGR;

    // W-component is only component in the liquid phase
    xmWL = 1.;

    rhoLR = liquid_phase.property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(variables, pos, t, dt);

    muLR = liquid_phase.property(MaterialPropertyLib::PropertyType::viscosity)
               .template value<double>(variables, pos, t, dt);

    lambdaLR =
        liquid_phase
            .property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .template value<double>(variables, pos, t, dt);

    rhoWLR = rhoLR;

    // specifi heat capacities
    auto const cpG =
        gas_phase
            .property(MaterialPropertyLib::PropertyType::specific_heat_capacity)
            .template value<double>(variables, pos, t, dt);

    auto const cpL =
        liquid_phase
            .property(MaterialPropertyLib::PropertyType::specific_heat_capacity)
            .template value<double>(variables, pos, t, dt);

    // specific phase enthalpies
    hG = cpG * T;
    hL = cpL * T;

    // specific inner energies
    uG = hG - pGR / rhoGR;
    uL = hL;
}
}  // namespace TH2M
}  // namespace ProcessLib
