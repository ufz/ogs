/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConstitutiveSetting.h"

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/MPL/Utils/FormKelvinVectorFromThermalExpansivity.h"
#include "MaterialLib/MPL/Utils/GetLiquidThermalExpansivity.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
void ConstitutiveSetting<DisplacementDim>::eval(
    double const t, double const dt,
    ParameterLib::SpatialPosition const& x_position,
    MaterialPropertyLib::Medium& medium, double const T, double const T_dot,
    GlobalDimVectorType const& grad_T, double const p_cap,
    double const p_cap_dot, GlobalDimVectorType const& grad_p_cap,
    KelvinVector const& eps_arg, KelvinVector const& eps_prev_arg)

{
    namespace MPL = MaterialPropertyLib;
    MPL::VariableArray variables;
    variables[static_cast<int>(MPL::Variable::capillary_pressure)] = p_cap;
    variables[static_cast<int>(MPL::Variable::phase_pressure)] = -p_cap;
    variables[static_cast<int>(MPL::Variable::temperature)] = T;

    MPL::VariableArray variables_prev;

    auto const& liquid_phase = medium.phase("AqueousLiquid");
    auto const& solid_phase = medium.phase("Solid");

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    eqU.eps.noalias() = eps_arg;

    auto& rho_LR = liquid_density;
    auto& S_L = saturation;
    auto const S_L_prev = saturation_prev;

    auto const alpha =
        medium.property(MPL::PropertyType::biot_coefficient)
            .template value<double>(variables, x_position, t, dt);

    double const dT = T_dot * dt;
    double const T_prev = T - dT;

    auto const C_el =
        eqU.computeElasticTangentStiffness(t, x_position, dt, T_prev, T);

    auto const beta_SR =
        (1 - alpha) / eqU.solid_material_.getBulkModulus(t, x_position, &C_el);
    variables[static_cast<int>(MPL::Variable::grain_compressibility)] = beta_SR;

    rho_LR = liquid_phase.property(MPL::PropertyType::density)
                 .template value<double>(variables, x_position, t, dt);

    S_L = medium.property(MPL::PropertyType::saturation)
              .template value<double>(variables, x_position, t, dt);
    variables[static_cast<int>(MPL::Variable::liquid_saturation)] = S_L;
    variables_prev[static_cast<int>(MPL::Variable::liquid_saturation)] =
        S_L_prev;

    auto const chi = [&medium, x_position, t, dt](double const S_L)
    {
        MPL::VariableArray vs;
        vs[static_cast<int>(MPL::Variable::liquid_saturation)] = S_L;
        return medium.property(MPL::PropertyType::bishops_effective_stress)
            .template value<double>(vs, x_position, t, dt);
    };

    double const chi_S_L = chi(S_L);

    variables[static_cast<int>(MPL::Variable::effective_pore_pressure)] =
        -chi_S_L * p_cap;

    {
        double const chi_S_L_prev = chi(S_L_prev);
        variables_prev[static_cast<int>(
            MPL::Variable::effective_pore_pressure)] =
            -chi_S_L_prev * (p_cap - p_cap_dot * dt);
    }

    // Set volumetric strain rate for the general case without swelling.
    variables[static_cast<int>(MPL::Variable::volumetric_strain)]
        .emplace<double>(Invariants::trace(eqU.eps));
    // TODO why not eps_prev? TODO check consistency
    variables_prev[static_cast<int>(MPL::Variable::volumetric_strain)]
        .emplace<double>(Invariants::trace(eps_prev_arg));

    auto& phi = porosity;
    {  // Porosity update
        variables_prev[static_cast<int>(MPL::Variable::porosity)] =
            porosity_prev;
        phi = medium.property(MPL::PropertyType::porosity)
                  .template value<double>(variables, variables_prev, x_position,
                                          t, dt);
        variables[static_cast<int>(MPL::Variable::porosity)] = phi;
    }

    if (alpha < phi)
    {
        OGS_FATAL(
            "ThermoRichardsMechanics: Biot-coefficient {} is smaller than "
            "porosity {} in element/integration point {}/{}.",
            alpha, phi, *x_position.getElementID(),
            *x_position.getIntegrationPoint());
    }

    // Swelling and possibly volumetric strain rate update.
    if (solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate))
    {
        // If there is swelling, compute it. Update volumetric strain rate,
        // s.t. it corresponds to the mechanical part only.
        eqU.sigma_sw = eqU.sigma_sw_prev;

        auto const sigma_sw_dot =
            MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
                MPL::formEigenTensor<3>(
                    solid_phase[MPL::PropertyType::swelling_stress_rate].value(
                        variables, variables_prev, x_position, t, dt)));
        eqU.sigma_sw += sigma_sw_dot * dt;

        // !!! Misusing volumetric strain for mechanical volumetric
        // strain just to update the transport porosity !!!
        std::get<double>(
            variables[static_cast<int>(MPL::Variable::volumetric_strain)]) +=
            identity2.transpose() * C_el.inverse() * eqU.sigma_sw;
        std::get<double>(variables_prev[static_cast<int>(
            MPL::Variable::volumetric_strain)]) +=
            identity2.transpose() * C_el.inverse() * eqU.sigma_sw_prev;
    }

    if (medium.hasProperty(MPL::PropertyType::transport_porosity))
    {
        variables_prev[static_cast<int>(MPL::Variable::transport_porosity)] =
            transport_porosity_prev;

        transport_porosity =
            medium.property(MPL::PropertyType::transport_porosity)
                .template value<double>(variables, variables_prev, x_position,
                                        t, dt);
        variables[static_cast<int>(MPL::Variable::transport_porosity)] =
            transport_porosity;
    }
    else
    {
        variables[static_cast<int>(MPL::Variable::transport_porosity)] = phi;
    }

    //
    // displacement equation, displacement part
    //

    // Consider also anisotropic thermal expansion.
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const
        solid_linear_thermal_expansivity_vector =
            MPL::formKelvinVectorFromThermalExpansivity<DisplacementDim>(
                solid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_expansivity)
                    .value(variables, x_position, t, dt));

    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const
        dthermal_strain = solid_linear_thermal_expansivity_vector * dT;

    eqU.eps_m.noalias() =
        eqU.eps_m_prev + eqU.eps - eqU.eps_prev - dthermal_strain;

    if (solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate))
    {
        eqU.eps_m.noalias() +=
            C_el.inverse() * (eqU.sigma_sw - eqU.sigma_sw_prev);
    }

    variables[static_cast<int>(
                  MaterialPropertyLib::Variable::mechanical_strain)]
        .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
            eqU.eps_m);

    eqU.stiffness_tensor =
        eqU.updateConstitutiveRelation(variables, t, x_position, dt, T_prev);
    auto& C = eqU.stiffness_tensor;

    {
        double const p_FR = -chi_S_L * p_cap;
        // p_SR
        variables[static_cast<int>(MPL::Variable::solid_grain_pressure)] =
            p_FR - Invariants::trace(eqU.sigma_eff) / (3 * (1 - phi));
    }

    auto const rho_SR =
        solid_phase.property(MPL::PropertyType::density)
            .template value<double>(variables, x_position, t, dt);
    dry_density_solid = (1 - phi) * rho_SR;

    eqU.sigma_total = eqU.sigma_eff + alpha * chi_S_L * p_cap * identity2;

    auto const b = process_data_.specific_body_force;

    {
        double const rho = rho_SR * (1 - phi) + S_L * phi * rho_LR;
        eqU.volumetric_body_force = rho * b;
    }

    //
    // displacement equation, pressure part
    //
    double const dS_L_dp_cap =
        medium.property(MPL::PropertyType::saturation)
            .template dValue<double>(variables,
                                     MPL::Variable::capillary_pressure,
                                     x_position, t, dt);

    {
        auto const dchi_dS_L =
            medium.property(MPL::PropertyType::bishops_effective_stress)
                .template dValue<double>(variables,
                                         MPL::Variable::liquid_saturation,
                                         x_position, t, dt);

        eqU.J_up_X_BTI2N = -alpha * (chi_S_L + dchi_dS_L * p_cap * dS_L_dp_cap);
        eqU.J_up_HT_V_N = phi * rho_LR * dS_L_dp_cap * b;
    }

    if (solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate))
    {
        using DimMatrix = Eigen::Matrix<double, 3, 3>;
        auto const dsigma_sw_dS_L =
            MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
                solid_phase.property(MPL::PropertyType::swelling_stress_rate)
                    .template dValue<DimMatrix>(
                        variables, variables_prev,
                        MPL::Variable::liquid_saturation, x_position, t, dt));

        eqU.J_up_BT_K_N = -dsigma_sw_dS_L * dS_L_dp_cap;
    }
    else
    {
        eqU.J_up_BT_K_N.setZero();
    }

    eqU.J_uT_BT_K_N =
        -C * solid_linear_thermal_expansivity_vector;  // TODO thermal stress?
    eqP.M_pu_X_BTI2N = S_L * rho_LR * alpha;

    //
    // pressure equation, pressure part.
    //
    double const k_rel =
        medium.property(MPL::PropertyType::relative_permeability)
            .template value<double>(variables, x_position, t, dt);

    viscosity = liquid_phase.property(MPL::PropertyType::viscosity)
                    .template value<double>(variables, x_position, t, dt);

    // Set mechanical variables for the intrinsic permeability model
    // For stress dependent permeability.

    // For stress dependent permeability.
    using SymmetricTensor =
        KelvinVector;  // same data type, but different semantics
    variables[static_cast<int>(MPL::Variable::total_stress)]
        .emplace<SymmetricTensor>(
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                eqU.sigma_total));

    variables[static_cast<int>(
        MaterialPropertyLib::Variable::equivalent_plastic_strain)] =
        eqU.material_state_variables_->getEquivalentPlasticStrain();

    auto const K_intrinsic = MPL::formEigenTensor<DisplacementDim>(
        medium.property(MPL::PropertyType::permeability)
            .value(variables, x_position, t, dt));

    GlobalDimMatrixType const Ki_over_mu = K_intrinsic / viscosity;
    GlobalDimMatrixType const rho_Ki_over_mu = rho_LR * Ki_over_mu;

    eqP.K_pp_Laplace = k_rel * rho_Ki_over_mu;

    const double alphaB_minus_phi = alpha - phi;
    double const drho_LR_dp =
        liquid_phase.property(MPL::PropertyType::density)
            .template dValue<double>(variables, MPL::Variable::phase_pressure,
                                     x_position, t, dt);

    auto const beta_LR = drho_LR_dp / rho_LR;

    double const a0 = alphaB_minus_phi * beta_SR;
    double const specific_storage_a_p = S_L * (phi * beta_LR + S_L * a0);
    double const specific_storage_a_S = phi - p_cap * S_L * a0;

    // Note: d beta_LR/d p is omitted because it is a small value.
    double const dspecific_storage_a_p_dp_cap =
        dS_L_dp_cap * (phi * beta_LR + 2 * S_L * a0);
    double const dspecific_storage_a_S_dp_cap =
        -a0 * (S_L + p_cap * dS_L_dp_cap);

    // secant derivative from time discretization for storage
    // use tangent, if secant is not available
    double const DeltaS_L_Deltap_cap =
        (p_cap_dot == 0) ? dS_L_dp_cap : (S_L - S_L_prev) / (dt * p_cap_dot);

    double const o_storage_p_a_p = rho_LR * specific_storage_a_p;
    eqP.storage_p_a_S_X_NTN =
        -rho_LR * specific_storage_a_S * DeltaS_L_Deltap_cap;
    eqP.J_pp_X_NTN = p_cap_dot * rho_LR * dspecific_storage_a_p_dp_cap;
    eqP.storage_p_a_S_Jpp_X_NTN =
        -rho_LR *
        ((S_L - S_L_prev) * dspecific_storage_a_S_dp_cap +
         specific_storage_a_S * dS_L_dp_cap) /
        dt;
    eqP.J_pp_X_BTI2NT_u_dot_N = -rho_LR * dS_L_dp_cap * alpha;

    double const dk_rel_dS_L =
        medium.property(MPL::PropertyType::relative_permeability)
            .template dValue<double>(
                variables, MPL::Variable::liquid_saturation, x_position, t, dt);

    eqP.J_pp_dNT_V_N =
        rho_Ki_over_mu * dk_rel_dS_L * dS_L_dp_cap * (grad_p_cap + rho_LR * b);
    GlobalDimVectorType const o_rhs_Darcy_gravity =
        rho_LR * (eqP.K_pp_Laplace * b);

    //
    // pressure equation, temperature part, thermal expansion.
    //
    {
        double const fluid_volumetric_thermal_expansion =
            phi * MPL::getLiquidThermalExpansivity(liquid_phase, variables,
                                                   rho_LR, x_position, t, dt);

        const double eff_thermal_expansion =
            alphaB_minus_phi *
                Invariants::trace(solid_linear_thermal_expansivity_vector) +
            fluid_volumetric_thermal_expansion;

        eqP.M_pT_X_NTN = -S_L * rho_LR * eff_thermal_expansion;
    }

    //
    // temperature equation.
    //
    GlobalDimVectorType advective_heat_flux_contribution_to_K_liquid =
        GlobalDimVectorType::Zero(DisplacementDim);
    {
        auto const specific_heat_capacity_fluid =
            liquid_phase.property(MaterialPropertyLib::specific_heat_capacity)
                .template value<double>(variables, x_position, t, dt);

        auto const specific_heat_capacity_solid =
            solid_phase
                .property(
                    MaterialPropertyLib::PropertyType::specific_heat_capacity)
                .template value<double>(variables, x_position, t, dt);

        // TODO real vs. "non-real" values
        double const volumetric_heat_capacity_liquid =
            rho_LR * specific_heat_capacity_fluid;

        // NOTE: Gas phase is not included.
        double const volumetric_heat_capacity_liquid_and_solid =
            rho_SR * specific_heat_capacity_solid * (1 - phi) +
            S_L * phi * volumetric_heat_capacity_liquid;

        eqT.M_TT_X_NTN = volumetric_heat_capacity_liquid_and_solid;

        eqT.K_TT_Laplace =
            MaterialPropertyLib::formEigenTensor<DisplacementDim>(
                medium
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .value(variables, x_position, t, dt));

        // TODO changed formula, double check
        v_darcy = Ki_over_mu * (k_rel * (grad_p_cap + rho_LR * b));

        // Unit is J / m^2 / s / K. It's not a heat flux, but related.
        advective_heat_flux_contribution_to_K_liquid =
            volumetric_heat_capacity_liquid * v_darcy;

        //
        // temperature equation, pressure part
        //
        eqT.K_Tp_NT_V_dN = -volumetric_heat_capacity_liquid * k_rel *
                           (Ki_over_mu.transpose() * grad_T);
        eqT.K_Tp_X_NTN = -volumetric_heat_capacity_liquid *
                         v_darcy.dot(grad_T) / k_rel * dk_rel_dS_L *
                         dS_L_dp_cap;
    }

    double volumetric_heat_capacity_vapor = 0;
    double storage_coefficient_by_water_vapor = 0;
    eqP.J_pT_X_dNTdN = 0;
    eqP.K_pp_X_dNTdN = 0;  // set to a proper value in the true-branch below
    eqT.K_TT_X_dNTdN = 0;  // set to a proper value in the true-branch below
    eqT.K_Tp_X_dNTdN = 0;  // set to a proper value in the true-branch below
    eqT.M_Tp_X_NTN = 0;    // set to a proper value in the true-branch below
    GlobalDimVectorType vapor_velocity =
        GlobalDimVectorType::Zero(DisplacementDim);
    if (liquid_phase.hasProperty(MPL::PropertyType::vapour_diffusion) &&
        S_L < 1.0)
    {
        variables[static_cast<int>(MPL::Variable::density)] = rho_LR;

        double const rho_wv =
            liquid_phase.property(MaterialPropertyLib::vapour_density)
                .template value<double>(variables, x_position, t, dt);

        double const drho_wv_dT =
            liquid_phase.property(MaterialPropertyLib::vapour_density)
                .template dValue<double>(variables, MPL::Variable::temperature,
                                         x_position, t, dt);
        double const drho_wv_dp =
            liquid_phase.property(MaterialPropertyLib::vapour_density)
                .template dValue<double>(variables,
                                         MPL::Variable::phase_pressure,
                                         x_position, t, dt);
        auto const f_Tv =
            liquid_phase
                .property(
                    MPL::PropertyType::thermal_diffusion_enhancement_factor)
                .template value<double>(variables, x_position, t, dt);

        variables[static_cast<int>(MPL::Variable::porosity)] = phi;
        double const D_v =
            liquid_phase.property(MPL::PropertyType::vapour_diffusion)
                .template value<double>(variables, x_position, t, dt);

        eqP.J_pT_X_dNTdN = f_Tv * D_v * drho_wv_dT;
        eqP.K_pp_X_dNTdN = D_v * drho_wv_dp;

        MPL::Phase const* const gas_phase =
            medium.hasPhase("Gas") ? &medium.phase("Gas") : nullptr;
        if (gas_phase &&
            gas_phase->hasProperty(MPL::PropertyType::specific_heat_capacity))
        {
            // Vapour velocity
            vapor_velocity =
                -(eqP.J_pT_X_dNTdN * grad_T - eqP.K_pp_X_dNTdN * grad_p_cap) /
                rho_LR;
            double const specific_heat_capacity_vapor =
                gas_phase
                    ->property(MaterialPropertyLib::PropertyType::
                                   specific_heat_capacity)
                    .template value<double>(variables, x_position, t, dt);

            volumetric_heat_capacity_vapor =
                rho_wv * specific_heat_capacity_vapor;
            eqT.M_TT_X_NTN += volumetric_heat_capacity_vapor * (1 - S_L) * phi;
        }

        storage_coefficient_by_water_vapor =
            phi * (rho_wv * dS_L_dp_cap + (1 - S_L) * drho_wv_dp);

        eqP.M_pT_X_NTN += phi * (1 - S_L) * drho_wv_dT;

        //
        // Latent heat term
        //
        if (liquid_phase.hasProperty(MPL::PropertyType::latent_heat))
        {
            double const factor = phi * (1 - S_L) / rho_LR;
            // The volumetric latent heat of vaporization of liquid water
            double const L0 =
                liquid_phase.property(MPL::PropertyType::latent_heat)
                    .template value<double>(variables, x_position, t, dt) *
                rho_LR;

            double const drho_LR_dT =
                liquid_phase.property(MPL::PropertyType::density)
                    .template dValue<double>(variables,
                                             MPL::Variable::temperature,
                                             x_position, t, dt);

            double const rho_wv_over_rho_L = rho_wv / rho_LR;

            // TODO += is not good
            eqT.M_TT_X_NTN +=
                factor * L0 * (drho_wv_dT - rho_wv_over_rho_L * drho_LR_dT);
            eqT.M_Tp_X_NTN =
                (factor * L0 * (drho_wv_dp - rho_wv_over_rho_L * drho_LR_dp) +
                 L0 * phi * rho_wv_over_rho_L * dS_L_dp_cap);
            eqT.K_TT_X_dNTdN = L0 * eqP.J_pT_X_dNTdN / rho_LR;
            eqT.K_Tp_X_dNTdN = L0 * eqP.K_pp_X_dNTdN / rho_LR;
        }
    }

    eqT.K_TT_NT_V_dN = advective_heat_flux_contribution_to_K_liquid +
                       vapor_velocity * volumetric_heat_capacity_vapor;

    eqP.storage_p_a_p_X_NTN =
        o_storage_p_a_p + storage_coefficient_by_water_vapor;
    eqP.rhs_p_dNT_V = -o_rhs_Darcy_gravity + eqP.J_pT_X_dNTdN * grad_T;
}

template struct ConstitutiveSetting<2>;
template struct ConstitutiveSetting<3>;

}  // namespace ProcessLib::ThermoRichardsMechanics
