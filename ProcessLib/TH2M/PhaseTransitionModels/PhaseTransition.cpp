/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "PhaseTransition.h"

#include "MaterialLib/PhysicalConstant.h"

namespace ProcessLib
{
namespace TH2M
{
PhaseTransition::PhaseTransition(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
    : PhaseTransitionModel(media),
      n_components_gas_{numberOfComponents(media, "Gas")},
      // Identifies and initialises the indices of the various components via
      // their properties
      gas_phase_vapour_component_index_{findComponentIndex(
          media, "Gas", MaterialPropertyLib::PropertyType::vapour_pressure)},
      // Dry air component is complement of vapour component index
      gas_phase_dry_air_component_index_{1 - gas_phase_vapour_component_index_},
      // The solute is the component of the liquid phase which has been assigned
      // the property `henry_coefficient'.
      liquid_phase_solute_component_index_{findComponentIndex(
          media, "AqueousLiquid",
          MaterialPropertyLib::PropertyType::henry_coefficient)},
      // The solvent is again the bitwise complement of the solute component
      liquid_phase_solvent_component_index_{
          1 - liquid_phase_solute_component_index_}
{
    DBUG("Create PhaseTransition constitutive model.");

    if (n_components_gas_ != 2)
    {
        OGS_FATAL(
            "The current implementation of PhaseTransitionModelEvaporation "
            "requires the specification of exactly two components in the "
            "gas "
            "phase.");
    }

    // It is always the first (begin) medium that holds fluid phases.
    auto const medium = media.begin()->second;
    auto const& gas_phase = medium->phase("Gas");
    auto const& liquid_phase = medium->phase("AqueousLiquid");

    // check for minimum requirement definitions in media object
    std::array const required_vapour_component_properties = {
        MaterialPropertyLib::vapour_pressure, MaterialPropertyLib::molar_mass,
        MaterialPropertyLib::specific_heat_capacity,
        MaterialPropertyLib::diffusion,
        MaterialPropertyLib::specific_latent_heat};

    std::array const required_dry_air_component_properties = {
        MaterialPropertyLib::molar_mass,
        MaterialPropertyLib::specific_heat_capacity};

    std::array const required_solute_component_properties = {
        MaterialPropertyLib::specific_heat_capacity,
        MaterialPropertyLib::henry_coefficient, MaterialPropertyLib::diffusion,
        MaterialPropertyLib::specific_latent_heat};

    std::array const required_solvent_component_properties = {
        MaterialPropertyLib::specific_heat_capacity};

    std::array const required_gas_properties = {MaterialPropertyLib::density,
                                                MaterialPropertyLib::viscosity};

    std::array const required_liquid_properties = {
        MaterialPropertyLib::density, MaterialPropertyLib::viscosity};

    checkRequiredProperties(
        gas_phase.component(gas_phase_vapour_component_index_),
        required_vapour_component_properties);
    checkRequiredProperties(
        gas_phase.component(gas_phase_dry_air_component_index_),
        required_dry_air_component_properties);
    checkRequiredProperties(
        liquid_phase.component(liquid_phase_solute_component_index_),
        required_solute_component_properties);
    checkRequiredProperties(
        liquid_phase.component(liquid_phase_solvent_component_index_),
        required_solvent_component_properties);
    checkRequiredProperties(gas_phase, required_gas_properties);
    checkRequiredProperties(liquid_phase, required_liquid_properties);
}

PhaseTransitionModelVariables PhaseTransition::updateConstitutiveVariables(
    PhaseTransitionModelVariables const& phase_transition_model_variables,
    const MaterialPropertyLib::Medium* medium,
    MaterialPropertyLib::VariableArray variables,
    ParameterLib::SpatialPosition pos, double const t, const double dt) const
{
    // primary variables
    auto const pGR = variables.gas_phase_pressure;
    auto const pCap = variables.capillary_pressure;
    auto const T = variables.temperature;
    auto const pLR = pGR - pCap;

    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& gas_phase = medium->phase("Gas");

    constexpr double R = MaterialLib::PhysicalConstant::IdealGasConstant;

    // identify vapour and air components for convenient access
    auto const& vapour_component =
        gas_phase.component(gas_phase_vapour_component_index_);
    auto const& dry_air_component =
        gas_phase.component(gas_phase_dry_air_component_index_);

    // specific latent heat (of evaporation)
    const auto dh_evap =
        vapour_component
            .property(MaterialPropertyLib::PropertyType::specific_latent_heat)
            .template value<double>(variables, pos, t, dt);

    // molar mass of evaporating component
    auto const M_W =
        vapour_component.property(MaterialPropertyLib::PropertyType::molar_mass)
            .template value<double>(variables, pos, t, dt);

    // provide evaporation enthalpy and molar mass of the evaporating
    // component in the variable array
    variables.enthalpy_of_evaporation = dh_evap;
    variables.molar_mass = M_W;

    // vapour pressure over flat interface
    const auto p_vap_flat =
        vapour_component
            .property(MaterialPropertyLib::PropertyType::vapour_pressure)
            .template value<double>(variables, pos, t, dt);

    const auto dp_vap_flat_dT =
        vapour_component
            .property(MaterialPropertyLib::PropertyType::vapour_pressure)
            .template dValue<double>(variables,
                                     MaterialPropertyLib::Variable::temperature,
                                     pos, t, dt);

    // molar mass of dry air component
    auto const M_C =
        dry_air_component
            .property(MaterialPropertyLib::PropertyType::molar_mass)
            .template value<double>(variables, pos, t, dt);

    // copy previous state before modification.
    PhaseTransitionModelVariables cv = phase_transition_model_variables;

    // Water pressure is passed to the VariableArray in order to calculate
    // the water density.
    variables.liquid_phase_pressure = pLR;

    // Concentration is initially zero to calculate the density of the pure
    // water phase, which is needed for the Kelvin-Laplace equation.
    variables.concentration = 0.;
    const auto rhoLR_0 =
        liquid_phase.property(MaterialPropertyLib::PropertyType::density)
            .template value<double>(variables, pos, t, dt);

    // Kelvin-Laplace correction for menisci
    const double K = pCap > 0. ? std::exp(-pCap * M_W / rhoLR_0 / R / T) : 1.;
    const double dK_dT = pCap > 0. ? pCap * M_W / rhoLR_0 / R / T / T * K : 0;
    const double dK_dpCap =
        pCap > 0. ? -M_W / rhoLR_0 / R / T * K
                  : 0.;  // rhoLR_0 is treated as a constant here. However, the
                         // resulting errors are very small and can be ignored.

    // vapour pressure inside porespace (== water partial pressure in gas
    // phase)
    cv.pWGR = p_vap_flat * K;
    auto const dpWGR_dT = dp_vap_flat_dT * K + p_vap_flat * dK_dT;
    auto const dpWGR_dpCap = p_vap_flat * dK_dpCap;

    // gas phase molar fractions
    auto const xnWG_min =
        1.e-12;  // Magic number; prevents the mass fraction of a component
                 // from ever becoming zero (which would cause the partial
                 // density to disappear and thus one of the Laplace terms
                 // of the mass balance on the diagonal of the local element
                 // matrix to be zero). The value is simply made up, seems
                 // reasonable.
    cv.xnWG = std::clamp(cv.pWGR / pGR, xnWG_min, 1. - xnWG_min);
    const double xnCG = 1. - cv.xnWG;

    // gas phase molar fraction derivatives
    auto const dxnWG_dpGR = -cv.pWGR / pGR / pGR;
    auto const dxnWG_dpCap = dpWGR_dpCap / pGR;
    auto const dxnWG_dT = dpWGR_dT / pGR;

    // molar mass of the gas phase as a mixture of 'air' and vapour
    auto const MG = xnCG * M_C + cv.xnWG * M_W;
    variables.molar_mass = MG;

    // gas phase mass fractions
    cv.xmWG = cv.xnWG * M_W / MG;
    const auto xmCG = 1. - cv.xmWG;

    auto const dxn_dxm_conversion = M_W * M_C / MG / MG;
    // gas phase mass fraction derivatives
    cv.dxmWG_dpGR = dxnWG_dpGR * dxn_dxm_conversion;
    cv.dxmWG_dpCap = dxnWG_dpCap * dxn_dxm_conversion;
    cv.dxmWG_dT = dxnWG_dT * dxn_dxm_conversion;

    // density of overall gas phase
    cv.rhoGR = gas_phase.property(MaterialPropertyLib::PropertyType::density)
                   .template value<double>(variables, pos, t, dt);

    // derivatives of average molar mass of the gas phase
    auto const dMG = M_W - M_C;
    auto const dMG_dpGR = dxnWG_dpGR * dMG;

    variables.molar_mass_derivative = dMG_dpGR;

    // Derivatives of the density of the (composite gas phase) and the
    // partial densities of its components. The density of the mixture can
    // be obtained via the property 'IdealGasLawBinaryMixture', for this
    // purpose the derivatives of the mean molar mass are passed in each
    // case via the variable_array.
    cv.drho_GR_dp_GR =
        gas_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(
                variables, MaterialPropertyLib::Variable::gas_phase_pressure,
                pos, t, dt);

    auto const dMG_dpCap = dxnWG_dpCap * dMG;
    variables.molar_mass_derivative = dMG_dpCap;
    cv.drho_GR_dp_cap =
        gas_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(
                variables, MaterialPropertyLib::Variable::capillary_pressure,
                pos, t, dt);

    auto const dMG_dT = dxnWG_dT * dMG;
    variables.molar_mass_derivative = dMG_dT;
    cv.drho_GR_dT =
        gas_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(variables,
                                     MaterialPropertyLib::Variable::temperature,
                                     pos, t, dt);

    // The derivatives of the partial densities of the gas phase are
    // hard-coded (they should remain so, as they are a fundamental part of
    // this evaporation model). By outsourcing the derivatives of the phase
    // density to the MPL, a constant phase density can also be assumed, the
    // derivatives of the partial densities are then unaffected and the
    // model is still consistent.
    cv.rhoCGR = xmCG * cv.rhoGR;
    cv.rhoWGR = cv.xmWG * cv.rhoGR;

    // 'Air'-component partial density derivatives
    cv.drho_C_GR_dp_GR = xmCG * cv.drho_GR_dp_GR - cv.dxmWG_dpGR * cv.rhoGR;
    cv.drho_C_GR_dp_cap = xmCG * cv.drho_GR_dp_cap - cv.dxmWG_dpCap * cv.rhoGR;
    cv.drho_C_GR_dT = xmCG * cv.drho_GR_dT - cv.dxmWG_dT * cv.rhoGR;

    // Vapour-component partial density derivatives
    cv.drho_W_GR_dp_GR = cv.xmWG * cv.drho_GR_dp_GR + cv.dxmWG_dpGR * cv.rhoGR;
    cv.drho_W_GR_dp_cap =
        cv.xmWG * cv.drho_GR_dp_cap + cv.dxmWG_dpCap * cv.rhoGR;
    cv.drho_W_GR_dT = cv.xmWG * cv.drho_GR_dT + cv.dxmWG_dT * cv.rhoGR;

    // specific heat capacities of dry air and vapour
    auto const cpCG =
        dry_air_component
            .property(MaterialPropertyLib::PropertyType::specific_heat_capacity)
            .template value<double>(variables, pos, t, dt);
    auto const cpWG =
        vapour_component
            .property(MaterialPropertyLib::PropertyType::specific_heat_capacity)
            .template value<double>(variables, pos, t, dt);

    // specific enthalpy of dry air and vapour components
    cv.hCG = cpCG * T;
    cv.hWG = cpWG * T + dh_evap;

    // specific enthalpy of gas phase
    cv.hG = xmCG * cv.hCG + cv.xmWG * cv.hWG;

    // specific inner energies of gas phase
    cv.uG = cv.hG - pGR / cv.rhoGR;

    // diffusion
    auto const tortuosity =
        medium->property(MaterialPropertyLib::PropertyType::tortuosity)
            .template value<double>(variables, pos, t, dt);

    auto const D_W_G_m =
        vapour_component.property(MaterialPropertyLib::PropertyType::diffusion)
            .template value<double>(variables, pos, t, dt);
    cv.diffusion_coefficient_vapour =
        tortuosity * D_W_G_m;  // Note here that D_W_G = D_C_G.

    variables.molar_fraction = xnCG;

    // gas phase viscosity
    cv.muGR = gas_phase.property(MaterialPropertyLib::PropertyType::viscosity)
                  .template value<double>(variables, pos, t, dt);

    // Dissolution part -- Liquid phase properties
    // -------------------------------------------

    // Reference to the gas component dissolved in the liquid phase
    auto const& solute_component =
        liquid_phase.component(liquid_phase_solute_component_index_);

    // The amount of dissolved gas is described by Henry's law. If no
    // dissolution is intended, the user must define the Henry coefficient
    // as 'constant 0'.
    //
    // Henry-Coefficient and derivatives
    auto const H =
        solute_component
            .property(MaterialPropertyLib::PropertyType::henry_coefficient)
            .template value<double>(variables, pos, t, dt);

    auto const dH_dT =
        solute_component
            .property(MaterialPropertyLib::PropertyType::henry_coefficient)
            .template dValue<double>(variables,
                                     MaterialPropertyLib::Variable::temperature,
                                     pos, t, dt);
    auto const dH_dpGR =
        solute_component
            .property(MaterialPropertyLib::PropertyType::henry_coefficient)
            .template dValue<double>(
                variables, MaterialPropertyLib::Variable::gas_phase_pressure,
                pos, t, dt);

    // Concentration of the dissolved gas as amount of substance of the
    // mixture component C related to the total volume of the liquid phase.
    auto const cCL = H * xnCG * pGR;
    // Fortunately for the developer, the signs of the derivatives of the
    // composition of binary mixtures are often opposed.
    auto const dxnCG_dpGR = -dxnWG_dpGR;
    auto const dxnCG_dT = -dxnWG_dT;

    auto const dcCL_dpGR = (dH_dpGR * xnCG + H * dxnCG_dpGR) * pGR + H * xnCG;
    auto const dcCL_dT = pGR * (dH_dT * xnCG + H * dxnCG_dT);

    // Density of pure liquid phase
    cv.rhoWLR = rhoLR_0;

    variables.concentration = cCL;
    // Liquid density including dissolved gas components. Attention! This
    // only works if the concentration of the C-component is taken into
    // account in the selected equation of state, e.g. via a (multi)-linear
    // equation of state. Should a constant density be assumed, no
    // dissolution will take place, since the composition of the water phase
    // is determined via the mass fractions (as the ratio of the densities
    // of solvent and solution). NB! This problem did not occur with the gas
    // phase because the composition there was determined via the molar
    // fractions (ratio of the partial pressures).
    cv.rhoLR = liquid_phase.property(MaterialPropertyLib::PropertyType::density)
                   .template value<double>(variables, pos, t, dt);
    variables.density = cv.rhoLR;

    // Gas component partial density in liquid phase
    cv.rhoCLR = cv.rhoLR - cv.rhoWLR;

    // liquid phase composition (mass fraction)
    cv.xmWL = std::clamp(cv.rhoWLR / cv.rhoLR, 0., 1.);
    auto const xmCL = 1. - cv.xmWL;

    // Attention! Usually a multi-linear equation of state is used to
    // determine the density of the solution. This requires independent
    // variables, but in this case the concentration is not independent, but
    // is determined via the equilibrium state. The exact derivation of the
    // density according to e.g. the pressure pLR is thus:
    //
    //  drho_d_pGR = rho_ref * (beta_pLR + betaC * dC_dpLR)
    //
    // instead of
    //
    //  d_rho_d_pGR =
    //     liquid_phase.property(MaterialPropertyLib::PropertyType::density)
    //         .template dValue<double>(
    //             variables,
    //             MaterialPropertyLib::Variable::liquid_phase_pressure,
    //             pos, t, dt);

    auto const rho_ref_betaP =
        liquid_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(
                variables, MaterialPropertyLib::Variable::liquid_phase_pressure,
                pos, t, dt);

    auto const rho_ref_betaT =
        liquid_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(variables,
                                     MaterialPropertyLib::Variable::temperature,
                                     pos, t, dt);

    auto const rho_ref_betaC =
        liquid_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(
                variables, MaterialPropertyLib::Variable::concentration, pos, t,
                dt);

    // liquid phase density derivatives
    auto const drhoLR_dpGR = rho_ref_betaP + rho_ref_betaC * dcCL_dpGR;
    auto const drhoLR_dpCap = -rho_ref_betaP;
    auto const drhoLR_dT = rho_ref_betaT + rho_ref_betaC * dcCL_dT;

    // solvent partial density derivatives
    auto const drhoWLR_dpGR = rho_ref_betaP;
    auto const drhoWLR_dpCap = -rho_ref_betaP;
    auto const drhoWLR_dT = rho_ref_betaT;

    // liquid phase mass fraction derivatives
    cv.dxmWL_dpGR = 1. / cv.rhoLR * (drhoWLR_dpGR - cv.xmWL * drhoLR_dpGR);
    cv.dxmWL_dpCap = 1. / cv.rhoLR * (drhoWLR_dpCap - cv.xmWL * drhoLR_dpCap);
    cv.dxmWL_dT = 1. / cv.rhoLR * (drhoWLR_dT - cv.xmWL * drhoLR_dT);

    // liquid phase molar fractions and derivatives
    cv.xnWL = cv.xmWL * M_C / (cv.xmWL * M_C + xmCL * M_W);

    // Reference to the pure liquid component
    auto const& solvent_component =
        liquid_phase.component(liquid_phase_solvent_component_index_);

    // specific heat capacities of liquid phase components
    auto const cpCL =
        solute_component
            .property(MaterialPropertyLib::PropertyType::specific_heat_capacity)
            .template value<double>(variables, pos, t, dt);
    auto const cpWL =
        solvent_component
            .property(MaterialPropertyLib::PropertyType::specific_heat_capacity)
            .template value<double>(variables, pos, t, dt);

    // specific heat of solution
    const auto dh_sol =
        solute_component
            .property(MaterialPropertyLib::PropertyType::specific_latent_heat)
            .template value<double>(variables, pos, t, dt);

    // specific enthalpy of liquid phase and its components
    cv.hCL = cpCL * T + dh_sol;
    cv.hWL = cpWL * T;
    cv.hL = xmCL * cv.hCL + cv.xmWL * cv.hWL;

    // specific inner energies of liquid phase
    cv.uL = cv.hL;

    // diffusion
    auto const D_C_L_m =
        solute_component.property(MaterialPropertyLib::PropertyType::diffusion)
            .template value<double>(variables, pos, t, dt);
    cv.diffusion_coefficient_solute =
        tortuosity * D_C_L_m;  // Note here that D_C_L = D_W_L.

    // liquid phase viscosity
    cv.muLR =
        liquid_phase.property(MaterialPropertyLib::PropertyType::viscosity)
            .template value<double>(variables, pos, t, dt);

    return cv;
}

}  // namespace TH2M
}  // namespace ProcessLib
