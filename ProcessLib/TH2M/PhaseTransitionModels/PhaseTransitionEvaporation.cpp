/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "PhaseTransitionEvaporation.h"

#include "MaterialLib/PhysicalConstant.h"

namespace ProcessLib
{
namespace TH2M
{
PhaseTransitionEvaporation::PhaseTransitionEvaporation(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
    : PhaseTransitionModel(media),
      n_components_gas_{numberOfComponents(media, "Gas")},
      gas_phase_vapour_component_index_{findComponentIndex(
          media, "Gas", MaterialPropertyLib::PropertyType::vapour_pressure)},
      // dry air component is complement of vapour component index
      gas_phase_dry_air_component_index_{gas_phase_vapour_component_index_ ^ 1}
{
    DBUG("Create PhaseTransitionEvaporation constitutive model.");

    if (n_components_gas_ != 2)
    {
        OGS_FATAL(
            "The current implementation of PhaseTransitionModelEvaporation "
            "requires the specification of exactly two components in the gas "
            "phase.");
    }

    // It is always the first (begin) medium that holds fluid phases.
    auto const medium = media.begin()->second;
    auto const& gas_phase = medium->phase("Gas");
    auto const& liquid_phase = medium->phase("AqueousLiquid");

    // check for minimum requirement definitions in media object
    std::array const required_vapour_component_properties = {
        MaterialPropertyLib::specific_latent_heat,
        MaterialPropertyLib::vapour_pressure, MaterialPropertyLib::molar_mass,
        MaterialPropertyLib::specific_heat_capacity,
        MaterialPropertyLib::diffusion};
    std::array const required_dry_air_component_properties = {
        MaterialPropertyLib::molar_mass,
        MaterialPropertyLib::specific_heat_capacity};
    std::array const required_liquid_properties = {
        MaterialPropertyLib::specific_heat_capacity};

    checkRequiredProperties(
        gas_phase.component(gas_phase_vapour_component_index_),
        required_vapour_component_properties);
    checkRequiredProperties(
        gas_phase.component(gas_phase_dry_air_component_index_),
        required_dry_air_component_properties);
    checkRequiredProperties(liquid_phase, required_liquid_properties);
}

PhaseTransitionModelVariables
PhaseTransitionEvaporation::updateConstitutiveVariables(
    PhaseTransitionModelVariables const& phase_transition_model_variables,
    const MaterialPropertyLib::Medium* medium,
    MaterialPropertyLib::VariableArray variables,
    ParameterLib::SpatialPosition pos, double const t, const double dt) const
{
    // primary variables
    auto const pGR = std::get<double>(variables[static_cast<int>(
        MaterialPropertyLib::Variable::phase_pressure)]);
    auto const pCap = std::get<double>(variables[static_cast<int>(
        MaterialPropertyLib::Variable::capillary_pressure)]);
    auto const T = std::get<double>(variables[static_cast<int>(
        MaterialPropertyLib::Variable::temperature)]);
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

    // provide evaporation enthalpy and molar mass of the evaporating component
    // in the variable array
    variables[static_cast<int>(
        MaterialPropertyLib::Variable::enthalpy_of_evaporation)] = dh_evap;
    variables[static_cast<int>(MaterialPropertyLib::Variable::molar_mass)] =
        M_W;

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

    variables[static_cast<int>(
        MaterialPropertyLib::Variable::liquid_phase_pressure)] = pLR;

    cv.rhoLR = liquid_phase.property(MaterialPropertyLib::PropertyType::density)
                   .template value<double>(variables, pos, t, dt);
    cv.rhoWLR = cv.rhoLR;

    // Kelvin-Laplace correction for menisci
    const double K = std::exp(-pCap * M_W / cv.rhoLR / R / T);
    const double dK_dT = pCap * M_W / cv.rhoLR / R / T / T * K;
    const double dK_dpCap =
        -M_W / cv.rhoLR / R / T *
        K;  // rhoLR is treated as a constant here, as well as the molar mass of
            // the liquid phase. However, the resulting errors are very small
            // and can be ignored.

    // vapour pressure inside porespace (== water partial pressure in gas phase)
    cv.pWGR = p_vap_flat * K;
    auto const dpWGR_dT = dp_vap_flat_dT * K + p_vap_flat * dK_dT;
    auto const dpWGR_dpCap = p_vap_flat * dK_dpCap;

    // gas phase molar fractions
    auto const xnWG_min =
        1.e-12;  // Magic number; prevents the mass fraction of a component from
                 // ever becoming zero (which would cause the partial density to
                 // disappear and thus one of the Laplace terms of the mass
                 // balance on the diagonal of the local element matrix to be
                 // zero). The value is simply made up, seems reasonable.
    cv.xnWG = std::clamp(cv.pWGR / pGR, xnWG_min, 1. - xnWG_min);
    const double xnCG = 1. - cv.xnWG;

    // gas phase molar fraction derivatives
    auto const dxnWG_dpGR = -cv.pWGR / pGR / pGR;
    auto const dxnWG_dpCap = dpWGR_dpCap / pGR;
    auto const dxnWG_dT = dpWGR_dT / pGR;

    // molar mass of the gas phase as a mixture of 'air' and vapour
    auto const MG = xnCG * M_C + cv.xnWG * M_W;
    variables[static_cast<int>(MaterialPropertyLib::Variable::molar_mass)] = MG;

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

    variables[static_cast<int>(
        MaterialPropertyLib::Variable::molar_mass_derivative)] = dMG_dpGR;

    // Derivatives of the density of the (composite gas phase) and the partial
    // densities of its components. The density of the mixture can be obtained
    // via the property 'IdealGasLawBinaryMixture', for this purpose the
    // derivatives of the mean molar mass are passed in each case via the
    // variable_array.
    cv.drho_GR_dp_GR =
        gas_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(
                variables, MaterialPropertyLib::Variable::phase_pressure, pos,
                t, dt);

    auto const dMG_dpCap = dxnWG_dpCap * dMG;
    variables[static_cast<int>(
        MaterialPropertyLib::Variable::molar_mass_derivative)] = dMG_dpCap;
    cv.drho_GR_dp_cap =
        gas_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(
                variables, MaterialPropertyLib::Variable::capillary_pressure,
                pos, t, dt);

    auto const dMG_dT = dxnWG_dT * dMG;
    variables[static_cast<int>(
        MaterialPropertyLib::Variable::molar_mass_derivative)] = dMG_dT;
    cv.drho_GR_dT =
        gas_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(variables,
                                     MaterialPropertyLib::Variable::temperature,
                                     pos, t, dt);

    // The derivatives of the partial densities of the gas phase are hard-coded
    // (they should remain so, as they are a fundamental part of this
    // evaporation model). By outsourcing the derivatives of the phase density
    // to the MPL, a constant phase density can also be assumed, the derivatives
    // of the partial densities are then unaffected and the model is still
    // consistent.
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

    // liquid phase composition
    cv.xmWL = 1.;
    cv.xnWL = 1.;

    // specific heat capacities of liquid phase
    auto const cpL =
        liquid_phase
            .property(MaterialPropertyLib::PropertyType::specific_heat_capacity)
            .template value<double>(variables, pos, t, dt);

    // specific enthalpy of liquid phase
    cv.hL = cpL * T;

    // specific inner energies of gas and liquid phases
    cv.uG = cv.hG - pGR / cv.rhoGR;
    cv.uL = cv.hL;

    // diffusion
    auto const tortuosity =
        1.0;  // Tortuosity formally belongs in the conversion from molecular
              // diffusion coefficient to effective diffusion coefficient. For
              // the moment, however, it is only determined by the coefficient
              // (i.e. by parameter 'diffusion' in the PRJ-file) itself.

    auto const D_W_G_m =
        vapour_component.property(MaterialPropertyLib::PropertyType::diffusion)
            .template value<double>(variables, pos, t, dt);
    cv.diffusion_coefficient_vapour =
        tortuosity * D_W_G_m;  // Note here that D_W_G = D_C_G.

    variables[static_cast<int>(MaterialPropertyLib::Variable::molar_fraction)] =
        xnCG;

    // gas phase viscosity
    cv.muGR = gas_phase.property(MaterialPropertyLib::PropertyType::viscosity)
                  .template value<double>(variables, pos, t, dt);

    // liquid phase viscosity
    cv.muLR =
        liquid_phase.property(MaterialPropertyLib::PropertyType::viscosity)
            .template value<double>(variables, pos, t, dt);

    return cv;
}

}  // namespace TH2M
}  // namespace ProcessLib
