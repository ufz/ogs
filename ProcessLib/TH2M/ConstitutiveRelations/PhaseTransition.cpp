// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "PhaseTransition.h"

#include "MaterialLib/PhysicalConstant.h"

namespace
{
int numberOfComponents(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media,
    std::string const& phase_name)
{
    // Always the first (begin) medium that holds fluid phases.
    auto const& medium = media.begin()->second;
    return medium->phase(phase_name).numberOfComponents();
}

int findComponentIndex(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media,
    std::string const& phase_name,
    MaterialPropertyLib::PropertyType property_type)
{
    // It is always the first (begin) medium that holds fluid phases.
    auto const& medium = media.begin()->second;
    auto const& phase = medium->phase(phase_name);

    // find the component for which the property 'property_type' is defined
    for (std::size_t c = 0; c < phase.numberOfComponents(); c++)
    {
        if (phase.component(c).hasProperty(property_type))
        {
            return c;
        }
    }

    // A lot of checks can (and should) be done to make sure that the right
    // components with the right properties are used. For example, the names
    // of the components can be compared to check that the name of the
    // evaporable component does not also correspond to the name of the
    // solute.

    OGS_FATAL(
        "PhaseTransitionModel: findComponentIndex() could not find the "
        "specified property type '{:s}' in phase '{:s}'.",
        MaterialPropertyLib::property_enum_to_string[property_type],
        phase_name);
}
}  // namespace

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
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

    // check for required medium properties
    std::array const required_medium_properties = {
        MaterialPropertyLib::tortuosity};
    checkRequiredProperties(*medium, required_medium_properties);

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

    std::array const required_gas_properties = {MaterialPropertyLib::density};

    std::array const required_liquid_properties = {
        MaterialPropertyLib::density};

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

void PhaseTransition::eval(SpaceTimeData const& x_t,
                           MediaData const& media_data,
                           GasPressureData const& p_GR,
                           CapillaryPressureData const& p_cap,
                           TemperatureData const& T_data,
                           PureLiquidDensityData const& rho_W_LR,
                           FluidEnthalpyData& fluid_enthalpy_data,
                           MassMoleFractionsData& mass_mole_fractions_data,
                           FluidDensityData& fluid_density_data,
                           VapourPartialPressureData& vapour_pressure_data,
                           ConstituentDensityData& constituent_density_data,
                           PhaseTransitionData& cv) const
{
    MaterialPropertyLib::VariableArray variables;

    // primary variables
    auto const pGR = p_GR.pG;
    auto const pCap = p_cap.pCap;
    auto const T = T_data.T;
    variables.gas_phase_pressure = pGR;
    variables.capillary_pressure = pCap;
    variables.temperature = T;

    auto const& liquid_phase = media_data.liquid;
    auto const& gas_phase = media_data.gas;

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
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    // molar mass of evaporating component
    auto const M_W =
        vapour_component.property(MaterialPropertyLib::PropertyType::molar_mass)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    // provide evaporation enthalpy and molar mass of the evaporating
    // component in the variable array
    variables.enthalpy_of_evaporation = dh_evap;
    variables.molar_mass = M_W;

    // vapour pressure over flat interface
    const auto p_vap_flat =
        vapour_component
            .property(MaterialPropertyLib::PropertyType::vapour_pressure)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    const auto dp_vap_flat_dT =
        vapour_component
            .property(MaterialPropertyLib::PropertyType::vapour_pressure)
            .template dValue<double>(variables,
                                     MaterialPropertyLib::Variable::temperature,
                                     x_t.x, x_t.t, x_t.dt);

    // molar mass of dry air component
    auto const M_C =
        dry_air_component
            .property(MaterialPropertyLib::PropertyType::molar_mass)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    // Water pressure is passed to the VariableArray in order to calculate
    // the water density.
    auto const pLR = pGR - pCap;
    variables.liquid_phase_pressure = pLR;

    // Kelvin-Laplace correction for menisci
    const double K =
        pCap > 0. ? std::exp(-pCap * M_W / rho_W_LR() / R / T) : 1.;
    const double dK_dT =
        pCap > 0. ? pCap * M_W / rho_W_LR() / R / T / T * K : 0;
    const double dK_dpCap =
        pCap > 0. ? -M_W / rho_W_LR() / R / T * K
                  : 0.;  // rho_W_LR is treated as a constant here. However, the
                         // resulting errors are very small and can be ignored.

    // vapour pressure inside porespace (== water partial pressure in gas
    // phase)
    vapour_pressure_data.pWGR = p_vap_flat * K;
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
    double const xnWG =
        std::clamp(vapour_pressure_data.pWGR / pGR, xnWG_min, 1. - xnWG_min);
    mass_mole_fractions_data.xnCG = 1. - xnWG;

    // gas phase molar fraction derivatives
    auto const dxnWG_dpGR = -vapour_pressure_data.pWGR / pGR / pGR;
    auto const dxnWG_dpCap = dpWGR_dpCap / pGR;
    auto const dxnWG_dT = dpWGR_dT / pGR;

    // molar mass of the gas phase as a mixture of 'air' and vapour
    auto const MG = mass_mole_fractions_data.xnCG * M_C + xnWG * M_W;
    variables.molar_mass = MG;

    // gas phase mass fractions
    double const xmWG = xnWG * M_W / MG;
    mass_mole_fractions_data.xmCG = 1. - xmWG;

    auto const dxn_dxm_conversion = M_W * M_C / MG / MG;
    // gas phase mass fraction derivatives
    cv.dxmWG_dpGR = dxnWG_dpGR * dxn_dxm_conversion;
    cv.dxmWG_dpCap = dxnWG_dpCap * dxn_dxm_conversion;
    cv.dxmWG_dT = dxnWG_dT * dxn_dxm_conversion;

    // density of overall gas phase
    fluid_density_data.rho_GR =
        gas_phase.property(MaterialPropertyLib::PropertyType::density)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

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
                x_t.x, x_t.t, x_t.dt);

    auto const dMG_dpCap = dxnWG_dpCap * dMG;
    variables.molar_mass_derivative = dMG_dpCap;
    cv.drho_GR_dp_cap =
        gas_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(
                variables, MaterialPropertyLib::Variable::capillary_pressure,
                x_t.x, x_t.t, x_t.dt);

    auto const dMG_dT = dxnWG_dT * dMG;
    variables.molar_mass_derivative = dMG_dT;
    cv.drho_GR_dT =
        gas_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(variables,
                                     MaterialPropertyLib::Variable::temperature,
                                     x_t.x, x_t.t, x_t.dt);

    // The derivatives of the partial densities of the gas phase are
    // hard-coded (they should remain so, as they are a fundamental part of
    // this evaporation model). By outsourcing the derivatives of the phase
    // density to the MPL, a constant phase density can also be assumed, the
    // derivatives of the partial densities are then unaffected and the
    // model is still consistent.
    constituent_density_data.rho_C_GR =
        mass_mole_fractions_data.xmCG * fluid_density_data.rho_GR;
    constituent_density_data.rho_W_GR = xmWG * fluid_density_data.rho_GR;

    // 'Air'-component partial density derivatives
    cv.drho_C_GR_dp_GR = mass_mole_fractions_data.xmCG * cv.drho_GR_dp_GR -
                         cv.dxmWG_dpGR * fluid_density_data.rho_GR;
    cv.drho_C_GR_dp_cap = mass_mole_fractions_data.xmCG * cv.drho_GR_dp_cap -
                          cv.dxmWG_dpCap * fluid_density_data.rho_GR;
    cv.drho_C_GR_dT = mass_mole_fractions_data.xmCG * cv.drho_GR_dT -
                      cv.dxmWG_dT * fluid_density_data.rho_GR;

    // Vapour-component partial density derivatives
    cv.drho_W_GR_dp_GR =
        xmWG * cv.drho_GR_dp_GR + cv.dxmWG_dpGR * fluid_density_data.rho_GR;
    cv.drho_W_GR_dp_cap =
        xmWG * cv.drho_GR_dp_cap + cv.dxmWG_dpCap * fluid_density_data.rho_GR;
    cv.drho_W_GR_dT =
        xmWG * cv.drho_GR_dT + cv.dxmWG_dT * fluid_density_data.rho_GR;

    // specific heat capacities of dry air and vapour
    auto const cpCG =
        dry_air_component
            .property(MaterialPropertyLib::PropertyType::specific_heat_capacity)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);
    auto const cpWG =
        vapour_component
            .property(MaterialPropertyLib::PropertyType::specific_heat_capacity)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    // specific enthalpy of dry air and vapour components
    cv.hCG = cpCG * T;
    cv.hWG = cpWG * T + dh_evap;

    // specific enthalpy of gas phase
    fluid_enthalpy_data.h_G =
        mass_mole_fractions_data.xmCG * cv.hCG + xmWG * cv.hWG;
    cv.dh_G_dT = 0;

    // specific inner energies of gas phase
    cv.uG = fluid_enthalpy_data.h_G - pGR / fluid_density_data.rho_GR;
    cv.du_G_dT = 0;
    cv.du_G_dp_GR = 0;

    // diffusion
    assert(media_data.tortuosity_prop != nullptr);
    auto const tortuosity = media_data.tortuosity_prop->template value<double>(
        variables, x_t.x, x_t.t, x_t.dt);

    auto const D_W_G_m =
        vapour_component.property(MaterialPropertyLib::PropertyType::diffusion)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);
    cv.diffusion_coefficient_vapour =
        tortuosity * D_W_G_m;  // Note here that D_W_G = D_C_G.

    variables.molar_fraction = mass_mole_fractions_data.xnCG;

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
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    auto const dH_dT =
        solute_component
            .property(MaterialPropertyLib::PropertyType::henry_coefficient)
            .template dValue<double>(variables,
                                     MaterialPropertyLib::Variable::temperature,
                                     x_t.x, x_t.t, x_t.dt);
    auto const dH_dpGR =
        solute_component
            .property(MaterialPropertyLib::PropertyType::henry_coefficient)
            .template dValue<double>(
                variables, MaterialPropertyLib::Variable::gas_phase_pressure,
                x_t.x, x_t.t, x_t.dt);

    // Concentration of the dissolved gas as amount of substance of the
    // mixture component C related to the total volume of the liquid phase.
    auto const cCL = H * mass_mole_fractions_data.xnCG * pGR;
    // Fortunately for the developer, the signs of the derivatives of the
    // composition of binary mixtures are often opposed.
    auto const dxnCG_dpGR = -dxnWG_dpGR;
    auto const dxnCG_dT = -dxnWG_dT;

    auto const dcCL_dpGR =
        (dH_dpGR * mass_mole_fractions_data.xnCG + H * dxnCG_dpGR) * pGR +
        H * mass_mole_fractions_data.xnCG;
    auto const dcCL_dT =
        pGR * (dH_dT * mass_mole_fractions_data.xnCG + H * dxnCG_dT);

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
    fluid_density_data.rho_LR =
        liquid_phase.property(MaterialPropertyLib::PropertyType::density)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);
    variables.density = fluid_density_data.rho_LR;

    // Gas component partial density in liquid phase
    constituent_density_data.rho_C_LR = fluid_density_data.rho_LR - rho_W_LR();

    // liquid phase composition (mass fraction)
    mass_mole_fractions_data.xmWL =
        std::clamp(rho_W_LR() / fluid_density_data.rho_LR, 0., 1.);
    auto const xmCL = 1. - mass_mole_fractions_data.xmWL;

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
    //             x_t.x, x_t.t, x_t.dt);

    auto const rho_ref_betaP =
        liquid_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(
                variables, MaterialPropertyLib::Variable::liquid_phase_pressure,
                x_t.x, x_t.t, x_t.dt);

    auto const rho_ref_betaT =
        liquid_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(variables,
                                     MaterialPropertyLib::Variable::temperature,
                                     x_t.x, x_t.t, x_t.dt);

    auto const rho_ref_betaC =
        liquid_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(
                variables, MaterialPropertyLib::Variable::concentration, x_t.x,
                x_t.t, x_t.dt);

    // liquid phase density derivatives
    auto const drhoLR_dpGR = rho_ref_betaP + rho_ref_betaC * dcCL_dpGR;
    auto const drhoLR_dpCap = -rho_ref_betaP;
    cv.drho_LR_dT = rho_ref_betaT + rho_ref_betaC * dcCL_dT;

    // solvent partial density derivatives
    auto const drhoWLR_dpGR = rho_ref_betaP;
    auto const drhoWLR_dpCap = -rho_ref_betaP;
    auto const drhoWLR_dT = rho_ref_betaT;

    // liquid phase mass fraction derivatives
    cv.dxmWL_dpGR =
        1. / fluid_density_data.rho_LR *
        (drhoWLR_dpGR - mass_mole_fractions_data.xmWL * drhoLR_dpGR);
    cv.dxmWL_dpCap =
        1. / fluid_density_data.rho_LR *
        (drhoWLR_dpCap - mass_mole_fractions_data.xmWL * drhoLR_dpCap);
    cv.dxmWL_dT = 1. / fluid_density_data.rho_LR *
                  (drhoWLR_dT - mass_mole_fractions_data.xmWL * cv.drho_LR_dT);

    // liquid phase molar fractions and derivatives
    mass_mole_fractions_data.xnWL =
        mass_mole_fractions_data.xmWL * M_C /
        (mass_mole_fractions_data.xmWL * M_C + xmCL * M_W);

    // Reference to the pure liquid component
    auto const& solvent_component =
        liquid_phase.component(liquid_phase_solvent_component_index_);

    // specific heat capacities of liquid phase components
    auto const cpCL =
        solute_component
            .property(MaterialPropertyLib::PropertyType::specific_heat_capacity)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);
    auto const cpWL =
        solvent_component
            .property(MaterialPropertyLib::PropertyType::specific_heat_capacity)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    // specific heat of solution
    const auto dh_sol =
        solute_component
            .property(MaterialPropertyLib::PropertyType::specific_latent_heat)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    // specific enthalpy of liquid phase and its components
    double const hCL = cpCL * T + dh_sol;
    double const hWL = cpWL * T;
    fluid_enthalpy_data.h_L = xmCL * hCL + mass_mole_fractions_data.xmWL * hWL;
    cv.dh_L_dT = 0;

    // specific inner energies of liquid phase
    cv.uL = fluid_enthalpy_data.h_L;
    cv.du_L_dT = 0;

    // diffusion
    auto const D_C_L_m =
        solute_component.property(MaterialPropertyLib::PropertyType::diffusion)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);
    cv.diffusion_coefficient_solute =
        tortuosity * D_C_L_m;  // Note here that D_C_L = D_W_L.

    // Some default initializations.
    cv.drho_LR_dp_LR = 0;
    cv.drho_W_LR_dp_GR = 0.;
    cv.drho_W_LR_dT = 0.;
    cv.drho_W_LR_dp_LR = 0.;
}

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
