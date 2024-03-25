/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/LU>

#include "ConstitutiveRelations/ConstitutiveData.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MathLib/EigenBlockMatrixView.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Reflection/ReflectionSetIPData.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "TH2MProcessData.h"

namespace ProcessLib
{
namespace TH2M
{
namespace MPL = MaterialPropertyLib;

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
TH2MLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                   DisplacementDim>::
    TH2MLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        TH2MProcessData<DisplacementDim>& process_data)
    : LocalAssemblerInterface<DisplacementDim>(
          e, integration_method, is_axially_symmetric, process_data)
{
    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();

    _ip_data.resize(n_integration_points);
    _secondary_data.N_u.resize(n_integration_points);

    auto const shape_matrices_u =
        NumLib::initShapeMatrices<ShapeFunctionDisplacement,
                                  ShapeMatricesTypeDisplacement,
                                  DisplacementDim>(e, is_axially_symmetric,
                                                   this->integration_method_);

    auto const shape_matrices_p =
        NumLib::initShapeMatrices<ShapeFunctionPressure,
                                  ShapeMatricesTypePressure, DisplacementDim>(
            e, is_axially_symmetric, this->integration_method_);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto& ip_data = _ip_data[ip];
        auto const& sm_u = shape_matrices_u[ip];
        ip_data.integration_weight =
            this->integration_method_.getWeightedPoint(ip).getWeight() *
            sm_u.integralMeasure * sm_u.detJ;

        ip_data.N_u = sm_u.N;
        ip_data.dNdx_u = sm_u.dNdx;

        ip_data.N_p = shape_matrices_p[ip].N;
        ip_data.dNdx_p = shape_matrices_p[ip].dNdx;

        _secondary_data.N_u[ip] = shape_matrices_u[ip].N;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::tuple<
    std::vector<ConstitutiveRelations::ConstitutiveData<DisplacementDim>>,
    std::vector<ConstitutiveRelations::ConstitutiveTempData<DisplacementDim>>>
TH2MLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                   DisplacementDim>::
    updateConstitutiveVariables(
        Eigen::VectorXd const& local_x, Eigen::VectorXd const& local_x_prev,
        double const t, double const dt,
        ConstitutiveRelations::ConstitutiveModels<DisplacementDim> const&
            models)
{
    [[maybe_unused]] auto const matrix_size =
        gas_pressure_size + capillary_pressure_size + temperature_size +
        displacement_size;

    assert(local_x.size() == matrix_size);

    auto const gas_pressure =
        local_x.template segment<gas_pressure_size>(gas_pressure_index);
    auto const capillary_pressure =
        local_x.template segment<capillary_pressure_size>(
            capillary_pressure_index);

    auto const temperature =
        local_x.template segment<temperature_size>(temperature_index);
    auto const temperature_prev =
        local_x_prev.template segment<temperature_size>(temperature_index);

    auto const displacement =
        local_x.template segment<displacement_size>(displacement_index);
    auto const displacement_prev =
        local_x_prev.template segment<displacement_size>(displacement_index);

    auto const& medium =
        *this->process_data_.media_map.getMedium(this->element_.getID());
    auto const& gas_phase = medium.phase("Gas");
    auto const& liquid_phase = medium.phase("AqueousLiquid");
    auto const& solid_phase = medium.phase("Solid");
    MediaData media_data{medium};

    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();

    std::vector<ConstitutiveRelations::ConstitutiveData<DisplacementDim>>
        ip_constitutive_data(n_integration_points);
    std::vector<ConstitutiveRelations::ConstitutiveTempData<DisplacementDim>>
        ip_constitutive_variables(n_integration_points);

    PhaseTransitionModelVariables ptmv;

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto& ip_data = _ip_data[ip];
        auto& ip_cv = ip_constitutive_variables[ip];
        auto& ip_cd = ip_constitutive_data[ip];
        auto& ip_out = this->output_data_[ip];
        auto& current_state = this->current_states_[ip];
        auto& prev_state = this->prev_states_[ip];

        auto const& Np = ip_data.N_p;
        auto const& NT = Np;
        auto const& Nu = ip_data.N_u;
        auto const& gradNu = ip_data.dNdx_u;
        auto const& gradNp = ip_data.dNdx_p;
        ParameterLib::SpatialPosition const pos{
            std::nullopt, this->element_.getID(), ip,
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    this->element_, Nu))};
        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                this->element_, Nu);

        double const T = NT.dot(temperature);
        double const T_prev = NT.dot(temperature_prev);
        TemperatureData const T_data{T, T_prev};
        double const pGR = Np.dot(gas_pressure);
        double const pCap = Np.dot(capillary_pressure);
        double const pLR = pGR - pCap;
        GlobalDimVectorType const gradpGR = gradNp * gas_pressure;
        GlobalDimVectorType const gradpCap = gradNp * capillary_pressure;
        GlobalDimVectorType const gradT = gradNp * temperature;

        // medium properties
        models.elastic_tangent_stiffness_model.eval({pos, t, dt}, T_data,
                                                    ip_cv.C_el_data);

        models.biot_model.eval({pos, t, dt}, media_data, ip_cv.biot_data);

        auto const Bu =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                gradNu, Nu, x_coord, this->is_axially_symmetric_);

        auto& eps = ip_out.eps_data.eps;
        eps.noalias() = Bu * displacement;
        models.S_L_model.eval({pos, t, dt}, media_data,
                              CapillaryPressureData{pCap},
                              current_state.S_L_data, ip_cv.dS_L_dp_cap);

        models.chi_S_L_model.eval({pos, t, dt}, media_data,
                                  current_state.S_L_data, ip_cv.chi_S_L);

        // solid phase compressibility
        models.beta_p_SR_model.eval({pos, t, dt}, ip_cv.biot_data,
                                    ip_cv.C_el_data, ip_cv.beta_p_SR);

        // If there is swelling stress rate, compute swelling stress.
        models.swelling_model.eval(
            {pos, t, dt}, media_data, ip_cv.C_el_data, current_state.S_L_data,
            prev_state.S_L_data, prev_state.swelling_data,
            current_state.swelling_data, ip_cv.swelling_data);

        // solid phase linear thermal expansion coefficient
        models.s_therm_exp_model.eval({pos, t, dt}, media_data,
                                      ip_cv.s_therm_exp_data);

        models.mechanical_strain_model.eval(
            T_data, ip_cv.s_therm_exp_data, ip_out.eps_data,
            Bu * displacement_prev, prev_state.mechanical_strain_data,
            ip_cv.swelling_data, current_state.mechanical_strain_data);

        models.s_mech_model.eval(
            {pos, t, dt}, T_data, current_state.mechanical_strain_data,
            prev_state.mechanical_strain_data, prev_state.eff_stress_data,
            current_state.eff_stress_data, this->material_states_[ip],
            ip_cd.s_mech_data, ip_cv.equivalent_plastic_strain_data);

        models.total_stress_model.eval(
            current_state.eff_stress_data, ip_cv.biot_data, ip_cv.chi_S_L,
            GasPressureData{pGR}, CapillaryPressureData{pCap},
            ip_cv.total_stress_data);

        models.permeability_model.eval(
            {pos, t, dt}, media_data, current_state.S_L_data,
            CapillaryPressureData{pCap}, T_data, ip_cv.total_stress_data,
            ip_out.eps_data, ip_cv.equivalent_plastic_strain_data,
            ip_out.permeability_data);

        MPL::VariableArray vars;
        MPL::VariableArray vars_prev;
        vars.temperature = T;
        vars.gas_phase_pressure = pGR;
        vars.capillary_pressure = pCap;
        vars.liquid_phase_pressure = pLR;

        // Set volumetric strain for the general case without swelling.
        vars.volumetric_strain = Invariants::trace(eps);

        vars.liquid_saturation = current_state.S_L_data.S_L;
        vars_prev.liquid_saturation = prev_state.S_L_data->S_L;

        auto const rho_ref_SR =
            solid_phase.property(MPL::PropertyType::density)
                .template value<double>(
                    vars, pos, t, std::numeric_limits<double>::quiet_NaN());

        double const T0 = this->process_data_.reference_temperature(t, pos)[0];
        double const delta_T(T - T0);
        ip_data.thermal_volume_strain =
            ip_cv.s_therm_exp_data.beta_T_SR * delta_T;

        // initial porosity
        auto const phi_0 = medium.property(MPL::PropertyType::porosity)
                               .template value<double>(vars, pos, t, dt);

        auto const phi_S_0 = 1. - phi_0;

#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
        auto const& m = Invariants::identity2;
        double const div_u = m.transpose() * eps;

        const double phi_S = phi_S_0 * (1. + ip_data.thermal_volume_strain -
                                        ip_cv.biot_data() * div_u);
#else   // NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
        const double phi_S = phi_S_0;
#endif  // NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION

        // porosity
        ip_data.phi = 1. - phi_S;
        vars.porosity = ip_data.phi;

        // solid phase density
#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
        auto const rhoSR = rho_ref_SR * (1. - ip_data.thermal_volume_strain +
                                         (ip_cv.biot_data() - 1.) * div_u);
#else   // NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
        auto const rhoSR = rho_ref_SR;
#endif  // NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION

        // constitutive model object as specified in process creation
        auto& ptm = *this->process_data_.phase_transition_model_;
        ptmv = ptm.updateConstitutiveVariables(ptmv, &medium, vars, pos, t, dt);
        auto const& c = ptmv;

        auto const phi_L = current_state.S_L_data.S_L * ip_data.phi;
        auto const phi_G = (1. - current_state.S_L_data.S_L) * ip_data.phi;

        // thermal conductivity
        ip_data.lambda = MaterialPropertyLib::formEigenTensor<DisplacementDim>(
            medium
                .property(
                    MaterialPropertyLib::PropertyType::thermal_conductivity)
                .value(vars, pos, t, dt));

        auto const cpS =
            solid_phase.property(MPL::PropertyType::specific_heat_capacity)
                .template value<double>(vars, pos, t, dt);
        ip_data.h_S = cpS * T;
        auto const u_S = ip_data.h_S;

        ip_data.rho_u_eff = phi_G * c.rhoGR * c.uG + phi_L * c.rhoLR * c.uL +
                            phi_S * rhoSR * u_S;

        ip_data.rho_G_h_G = phi_G * c.rhoGR * c.hG;
        ip_data.rho_L_h_L = phi_L * c.rhoLR * c.hL;
        ip_data.rho_S_h_S = phi_S * rhoSR * ip_data.h_S;

        ip_data.muGR = c.muGR;
        ip_data.muLR = c.muLR;

        ip_data.rhoGR = c.rhoGR;
        ip_data.rhoLR = c.rhoLR;
        ip_data.rhoSR = rhoSR;

        ip_data.rhoCGR = c.rhoCGR;
        ip_data.rhoCLR = c.rhoCLR;
        ip_data.rhoWGR = c.rhoWGR;
        ip_data.rhoWLR = c.rhoWLR;

        ip_data.dxmWG_dpGR = c.dxmWG_dpGR;
        ip_data.dxmWG_dpCap = c.dxmWG_dpCap;
        ip_data.dxmWG_dT = c.dxmWG_dT;

        ip_data.dxmWL_dpGR = c.dxmWL_dpGR;
        ip_data.dxmWL_dpCap = c.dxmWL_dpCap;
        ip_data.dxmWL_dT = c.dxmWL_dT;

        ip_data.dxmWL_dpLR = c.dxmWL_dpLR;

        // for variable output
        ip_data.xnCG = 1. - c.xnWG;
        ip_data.xmCG = 1. - c.xmWG;
        ip_data.xmWG = c.xmWG;
        ip_data.xmWL = c.xmWL;
        auto const xmCL = 1. - c.xmWL;

        ip_data.diffusion_coefficient_vapour = c.diffusion_coefficient_vapour;
        ip_data.diffusion_coefficient_solute = c.diffusion_coefficient_solute;

        ip_data.h_G = c.hG;
        ip_data.h_CG = c.hCG;
        ip_data.h_WG = c.hWG;
        ip_data.h_L = c.hL;
        ip_data.pWGR = c.pWGR;

        const GlobalDimVectorType gradxmWG = ip_data.dxmWG_dpGR * gradpGR +
                                             ip_data.dxmWG_dpCap * gradpCap +
                                             ip_data.dxmWG_dT * gradT;
        const GlobalDimVectorType gradxmCG = -gradxmWG;

        const GlobalDimVectorType gradxmWL = ip_data.dxmWL_dpGR * gradpGR +
                                             ip_data.dxmWL_dpCap * gradpCap +
                                             ip_data.dxmWL_dT * gradT;
        const GlobalDimVectorType gradxmCL = -gradxmWL;

        // Todo: factor -phiAlpha / xmZetaAlpha * DZetaAlpha can be evaluated in
        // the respective phase transition model, here only the multiplication
        // with the gradient of the mass fractions should take place.

        ip_data.d_CG = ip_data.xmCG == 0.
                           ? 0. * gradxmCG  // Keep d_CG's dimension and prevent
                                            // division by zero
                           : -phi_G / ip_data.xmCG *
                                 ip_data.diffusion_coefficient_vapour *
                                 gradxmCG;

        ip_data.d_WG = ip_data.xmWG == 0.
                           ? 0. * gradxmWG  // Keep d_WG's dimension and prevent
                                            // division by zero
                           : -phi_G / ip_data.xmWG *
                                 ip_data.diffusion_coefficient_vapour *
                                 gradxmWG;

        ip_data.d_CL = xmCL == 0. ? 0. * gradxmCL  // Keep d_CL's dimension and
                                                   // prevent division by zero
                                  : -phi_L / xmCL *
                                        ip_data.diffusion_coefficient_solute *
                                        gradxmCL;

        ip_data.d_WL = ip_data.xmWL == 0.
                           ? 0. * gradxmWL  // Keep d_WG's dimension and prevent
                                            // division by zero
                           : -phi_L / ip_data.xmWL *
                                 ip_data.diffusion_coefficient_solute *
                                 gradxmWL;

        // ---------------------------------------------------------------------
        // Derivatives for Jacobian
        // ---------------------------------------------------------------------
        auto const drho_LR_dT =
            liquid_phase.property(MPL::PropertyType::density)
                .template dValue<double>(vars, MPL::Variable::temperature, pos,
                                         t, dt);
        auto const drho_SR_dT =
            solid_phase.property(MPL::PropertyType::density)
                    .template dValue<double>(vars, MPL::Variable::temperature,
                                             pos, t, dt)
#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
                * (1. - ip_data.thermal_volume_strain +
                   (ip_cv.biot_data() - 1.) * div_u) -
            rho_ref_SR * ip_cv.s_therm_exp_data.beta_T_SR
#endif
            ;

        // porosity
        auto const dphi_0_dT =
            medium[MPL::PropertyType::porosity].template dValue<double>(
                vars, MPL::Variable::temperature, pos, t, dt);

        auto const dphi_S_0_dT = -dphi_0_dT;
        const double dphi_S_dT = dphi_S_0_dT
#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
                                     * (1. + ip_data.thermal_volume_strain -
                                        ip_cv.biot_data() * div_u) +
                                 phi_S_0 * ip_cv.s_therm_exp_data.beta_T_SR
#endif
            ;

        ip_cv.drho_u_eff_dT =
            phi_G * c.drho_GR_dT * c.uG + phi_G * c.rhoGR * c.du_G_dT +
            phi_L * drho_LR_dT * c.uL + phi_L * c.rhoLR * c.du_L_dT +
            phi_S * drho_SR_dT * u_S + phi_S * rhoSR * cpS +
            dphi_S_dT * rhoSR * u_S;

        // ds_L_dp_GR = 0;
        // ds_G_dp_GR = -ds_L_dp_GR;
        double const ds_G_dp_cap = -ip_cv.dS_L_dp_cap();

        // dphi_G_dp_GR = -ds_L_dp_GR * ip_data.phi = 0;
        double const dphi_G_dp_cap = -ip_cv.dS_L_dp_cap() * ip_data.phi;
        // dphi_L_dp_GR = ds_L_dp_GR * ip_data.phi = 0;
        double const dphi_L_dp_cap = ip_cv.dS_L_dp_cap() * ip_data.phi;

        auto const lambdaGR =
            gas_phase.hasProperty(MPL::PropertyType::thermal_conductivity)
                ? MPL::formEigenTensor<DisplacementDim>(
                      gas_phase
                          .property(MPL::PropertyType::thermal_conductivity)
                          .value(vars, pos, t, dt))
                : MPL::formEigenTensor<DisplacementDim>(0.);

        auto const dlambda_GR_dT =
            gas_phase.hasProperty(MPL::PropertyType::thermal_conductivity)
                ? MPL::formEigenTensor<DisplacementDim>(
                      gas_phase[MPL::PropertyType::thermal_conductivity].dValue(
                          vars, MPL::Variable::temperature, pos, t, dt))
                : MPL::formEigenTensor<DisplacementDim>(0.);

        auto const lambdaLR =
            liquid_phase.hasProperty(MPL::PropertyType::thermal_conductivity)
                ? MPL::formEigenTensor<DisplacementDim>(
                      liquid_phase
                          .property(MPL::PropertyType::thermal_conductivity)
                          .value(vars, pos, t, dt))
                : MPL::formEigenTensor<DisplacementDim>(0.);

        auto const dlambda_LR_dT =
            liquid_phase.hasProperty(MPL::PropertyType::thermal_conductivity)
                ? MPL::formEigenTensor<DisplacementDim>(
                      liquid_phase[MPL::PropertyType::thermal_conductivity]
                          .dValue(vars, MPL::Variable::temperature, pos, t, dt))
                : MPL::formEigenTensor<DisplacementDim>(0.);

        auto const lambdaSR =
            solid_phase.hasProperty(MPL::PropertyType::thermal_conductivity)
                ? MPL::formEigenTensor<DisplacementDim>(
                      solid_phase
                          .property(MPL::PropertyType::thermal_conductivity)
                          .value(vars, pos, t, dt))
                : MPL::formEigenTensor<DisplacementDim>(0.);

        auto const dlambda_SR_dT =
            solid_phase.hasProperty(MPL::PropertyType::thermal_conductivity)
                ? MPL::formEigenTensor<DisplacementDim>(
                      solid_phase[MPL::PropertyType::thermal_conductivity]
                          .dValue(vars, MPL::Variable::temperature, pos, t, dt))
                : MPL::formEigenTensor<DisplacementDim>(0.);

        ip_cv.dlambda_dp_cap =
            dphi_G_dp_cap * lambdaGR + dphi_L_dp_cap * lambdaLR;

        ip_cv.dlambda_dT = phi_G * dlambda_GR_dT + phi_L * dlambda_LR_dT +
                           phi_S * dlambda_SR_dT + dphi_S_dT * lambdaSR;

        // From p_LR = p_GR - p_cap it follows for
        // drho_LR/dp_GR = drho_LR/dp_LR * dp_LR/dp_GR
        //               = drho_LR/dp_LR * (dp_GR/dp_GR - dp_cap/dp_GR)
        //               = drho_LR/dp_LR * (1 - 0)
        double const drho_LR_dp_GR = c.drho_LR_dp_LR;
        double const drho_LR_dp_cap = -c.drho_LR_dp_LR;
        // drho_GR_dp_cap = 0;

        ip_cv.drho_h_eff_dp_GR =
            /*(dphi_G_dp_GR = 0) * c.rhoGR * c.hG +*/ phi_G * c.drho_GR_dp_GR *
                c.hG +
            /*(dphi_L_dp_GR = 0) * c.rhoLR * c.hL +*/ phi_L * drho_LR_dp_GR *
                c.hL;
        ip_cv.drho_h_eff_dp_cap = dphi_G_dp_cap * c.rhoGR * c.hG +
                                  /*phi_G * (drho_GR_dp_cap = 0) * c.hG +*/
                                  dphi_L_dp_cap * c.rhoLR * c.hL +
                                  phi_L * drho_LR_dp_cap * c.hL;

        // TODO (naumov) Extend for temperature dependent porosities.
        constexpr double dphi_G_dT = 0;
        constexpr double dphi_L_dT = 0;
        ip_cv.drho_h_eff_dT =
            dphi_G_dT * c.rhoGR * c.hG + phi_G * c.drho_GR_dT * c.hG +
            phi_G * c.rhoGR * c.dh_G_dT + dphi_L_dT * c.rhoLR * c.hL +
            phi_L * drho_LR_dT * c.hL + phi_L * c.rhoLR * c.dh_L_dT +
            dphi_S_dT * rhoSR * ip_data.h_S + phi_S * drho_SR_dT * ip_data.h_S +
            phi_S * rhoSR * cpS;

        ip_cv.drho_u_eff_dp_GR =
            /*(dphi_G_dp_GR = 0) * c.rhoGR * c.uG +*/
            phi_G * c.drho_GR_dp_GR * c.uG + phi_G * c.rhoGR * c.du_G_dp_GR +
            /*(dphi_L_dp_GR = 0) * c.rhoLR * c.uL +*/
            phi_L * drho_LR_dp_GR * c.uL + phi_L * c.rhoLR * c.du_L_dp_GR;

        ip_cv.drho_u_eff_dp_cap = dphi_G_dp_cap * c.rhoGR * c.uG +
                                  /*phi_G * (drho_GR_dp_cap = 0) * c.uG +*/
                                  dphi_L_dp_cap * c.rhoLR * c.uL +
                                  phi_L * drho_LR_dp_cap * c.uL +
                                  phi_L * c.rhoLR * c.du_L_dp_cap;

        auto const& b = this->process_data_.specific_body_force;
        GlobalDimMatrixType const k_over_mu_G =
            ip_out.permeability_data.Ki * ip_out.permeability_data.k_rel_G /
            ip_data.muGR;
        GlobalDimMatrixType const k_over_mu_L =
            ip_out.permeability_data.Ki * ip_out.permeability_data.k_rel_L /
            ip_data.muLR;

        // dk_over_mu_G_dp_GR = ip_out.permeability_data.Ki *
        //                      ip_out.permeability_data.dk_rel_G_dS_L *
        //                      (ds_L_dp_GR = 0) / ip_data.muGR = 0;
        // dk_over_mu_L_dp_GR = ip_out.permeability_data.Ki *
        //                      ip_out.permeability_data.dk_rel_L_dS_L *
        //                      (ds_L_dp_GR = 0) / ip_data.muLR = 0;
        ip_cv.dk_over_mu_G_dp_cap = ip_out.permeability_data.Ki *
                                    ip_out.permeability_data.dk_rel_G_dS_L *
                                    ip_cv.dS_L_dp_cap() / ip_data.muGR;
        ip_cv.dk_over_mu_L_dp_cap = ip_out.permeability_data.Ki *
                                    ip_out.permeability_data.dk_rel_L_dS_L *
                                    ip_cv.dS_L_dp_cap() / ip_data.muLR;

        ip_data.w_GS = k_over_mu_G * c.rhoGR * b - k_over_mu_G * gradpGR;
        ip_data.w_LS = k_over_mu_L * gradpCap + k_over_mu_L * c.rhoLR * b -
                       k_over_mu_L * gradpGR;

        ip_cv.drho_GR_h_w_eff_dp_GR_Npart =
            c.drho_GR_dp_GR * c.hG * ip_data.w_GS +
            c.rhoGR * c.hG * k_over_mu_G * c.drho_GR_dp_GR * b;
        ip_cv.drho_GR_h_w_eff_dp_GR_gradNpart =
            -c.rhoGR * c.hG * k_over_mu_G - c.rhoLR * c.hL * k_over_mu_L;

        ip_cv.drho_LR_h_w_eff_dp_cap_Npart =
            -drho_LR_dp_cap * c.hL * ip_data.w_LS -
            c.rhoLR * c.hL * k_over_mu_L * drho_LR_dp_cap * b;
        ip_cv.drho_LR_h_w_eff_dp_cap_gradNpart =
            // TODO (naumov) why the minus sign??????
            -c.rhoLR * c.hL * k_over_mu_L;

        ip_cv.drho_GR_h_w_eff_dT = c.drho_GR_dT * c.hG * ip_data.w_GS +
                                   c.rhoGR * c.dh_G_dT * ip_data.w_GS +
                                   drho_LR_dT * c.hL * ip_data.w_LS +
                                   c.rhoLR * c.dh_L_dT * ip_data.w_LS;
        // TODO (naumov) + k_over_mu_G * drho_GR_dT * b + k_over_mu_L *
        // drho_LR_dT * b

        // Derivatives of s_G * rho_C_GR_dot + s_L * rho_C_LR_dot abbreviated
        // here with S_rho_C_eff.
        double const s_L = current_state.S_L_data.S_L;
        double const s_G = 1. - s_L;
        double const rho_C_FR = s_G * ip_data.rhoCGR + s_L * ip_data.rhoCLR;
        double const rho_W_FR = s_G * ip_data.rhoWGR + s_L * ip_data.rhoWLR;
        // TODO (naumov) Extend for partially saturated media.
        constexpr double drho_C_GR_dp_cap = 0;
        if (dt == 0.)
        {
            ip_cv.dfC_3a_dp_GR = 0.;
            ip_cv.dfC_3a_dp_cap = 0.;
            ip_cv.dfC_3a_dT = 0.;
        }
        else
        {
            double const rho_C_GR_dot =
                (ip_data.rhoCGR - ip_data.rhoCGR_prev) / dt;
            double const rho_C_LR_dot =
                (ip_data.rhoCLR - ip_data.rhoCLR_prev) / dt;
            ip_cv.dfC_3a_dp_GR =
                /*(ds_G_dp_GR = 0) * rho_C_GR_dot +*/ s_G * c.drho_C_GR_dp_GR /
                    dt +
                /*(ds_L_dp_GR = 0) * rho_C_LR_dot +*/ s_L * c.drho_C_LR_dp_GR /
                    dt;
            ip_cv.dfC_3a_dp_cap = ds_G_dp_cap * rho_C_GR_dot +
                                  s_G * drho_C_GR_dp_cap / dt +
                                  ip_cv.dS_L_dp_cap() * rho_C_LR_dot -
                                  s_L * c.drho_C_LR_dp_LR / dt;
            ip_cv.dfC_3a_dT =
                s_G * c.drho_C_GR_dT / dt + s_L * c.drho_C_LR_dT / dt;
        }

        double const drho_C_FR_dp_GR =
            /*(ds_G_dp_GR = 0) * ip_data.rhoCGR +*/ s_G * c.drho_C_GR_dp_GR +
            /*(ds_L_dp_GR = 0) * ip_data.rhoCLR +*/ s_L * c.drho_C_LR_dp_GR;
        ip_cv.dfC_4_MCpG_dp_GR = drho_C_FR_dp_GR *
                                 (ip_cv.biot_data() - ip_data.phi) *
                                 ip_cv.beta_p_SR();

        double const drho_C_FR_dT = s_G * c.drho_C_GR_dT + s_L * c.drho_C_LR_dT;
        ip_cv.dfC_4_MCpG_dT =
            drho_C_FR_dT * (ip_cv.biot_data() - ip_data.phi) * ip_cv.beta_p_SR()
#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
            - rho_C_FR * ip_data.dphi_dT * ip_cv.beta_p_SR()
#endif
            ;

        ip_cv.dfC_4_MCT_dT = drho_C_FR_dT * (ip_cv.biot_data() - ip_data.phi) *
                                 ip_cv.s_therm_exp_data.beta_T_SR
#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
                             +
                             rho_C_FR * (ip_cv.biot_data() - ip_data.dphi_dT) *
                                 ip_cv.s_therm_exp_data.beta_T_SR
#endif
            ;

        ip_cv.dfC_4_MCu_dT = drho_C_FR_dT * ip_cv.biot_data();

        ip_cv.dfC_2a_dp_GR = -ip_data.phi * c.drho_C_GR_dp_GR -
                             drho_C_FR_dp_GR * pCap *
                                 (ip_cv.biot_data() - ip_data.phi) *
                                 ip_cv.beta_p_SR();

        double const drho_C_FR_dp_cap =
            ds_G_dp_cap * ip_data.rhoCGR + s_G * drho_C_GR_dp_cap +
            ip_cv.dS_L_dp_cap() * ip_data.rhoCLR - s_L * c.drho_C_LR_dp_LR;

        ip_cv.dfC_2a_dp_cap =
            ip_data.phi * (-c.drho_C_LR_dp_LR - drho_C_GR_dp_cap) -
            drho_C_FR_dp_cap * pCap * (ip_cv.biot_data() - ip_data.phi) *
                ip_cv.beta_p_SR() +
            rho_C_FR * (ip_cv.biot_data() - ip_data.phi) * ip_cv.beta_p_SR();

        ip_cv.dfC_2a_dT =
#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
            ip_data.dphi_dT * (ip_data.rhoCLR - ip_data.rhoCGR) +
#endif
            ip_data.phi * (c.drho_C_LR_dT - c.drho_C_GR_dT) -
            drho_C_FR_dT * pCap * (ip_cv.biot_data() - ip_data.phi) *
                ip_cv.beta_p_SR()
#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
            + rho_C_FR * pCap * ip_data.dphi_dT * ip_cv.beta_p_SR()
#endif
            ;

        ip_cv.dadvection_C_dp_GR = c.drho_C_GR_dp_GR * k_over_mu_G
                                   // + rhoCGR * (dk_over_mu_G_dp_GR = 0)
                                   // + rhoCLR * (dk_over_mu_L_dp_GR = 0)
                                   + c.drho_C_LR_dp_GR * k_over_mu_L;

        ip_cv.dadvection_C_dp_cap =
            //(drho_C_GR_dp_cap = 0) * k_over_mu_G
            ip_data.rhoCGR * ip_cv.dk_over_mu_G_dp_cap +
            (-c.drho_C_LR_dp_LR) * k_over_mu_L +
            ip_data.rhoCLR * ip_cv.dk_over_mu_L_dp_cap;

        ip_cv.dfC_4_LCpG_dT =
            c.drho_C_GR_dT * k_over_mu_G + c.drho_C_LR_dT * k_over_mu_L
            // + ip_cv.ddiffusion_C_p_dT TODO (naumov)
            ;

        double const drho_W_FR_dp_GR =
            /*(ds_G_dp_GR = 0) * ip_data.rhoWGR +*/ s_G * c.drho_W_GR_dp_GR +
            /*(ds_L_dp_GR = 0) * ip_data.rhoWLR +*/ s_L * c.drho_W_LR_dp_GR;
        double const drho_W_FR_dp_cap =
            ds_G_dp_cap * ip_data.rhoWGR + s_G * c.drho_W_GR_dp_cap +
            ip_cv.dS_L_dp_cap() * ip_data.rhoWLR - s_L * c.drho_W_LR_dp_LR;
        double const drho_W_FR_dT = s_G * c.drho_W_GR_dT + s_L * c.drho_W_LR_dT;

        ip_cv.dfW_2a_dp_GR =
            ip_data.phi * (c.drho_W_LR_dp_GR - c.drho_W_GR_dp_GR);
        ip_cv.dfW_2b_dp_GR = drho_W_FR_dp_GR * pCap *
                             (ip_cv.biot_data() - ip_data.phi) *
                             ip_cv.beta_p_SR();
        ip_cv.dfW_2a_dp_cap =
            ip_data.phi * (-c.drho_W_LR_dp_LR - c.drho_W_GR_dp_cap);
        ip_cv.dfW_2b_dp_cap =
            drho_W_FR_dp_cap * pCap * (ip_cv.biot_data() - ip_data.phi) *
                ip_cv.beta_p_SR() +
            rho_W_FR * (ip_cv.biot_data() - ip_data.phi) * ip_cv.beta_p_SR();

        ip_cv.dfW_2a_dT =
#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
            ip_data.dphi_dT * (ip_data.rhoWLR - ip_data.rhoWGR) +
#endif
            ip_data.phi * (c.drho_W_LR_dT - c.drho_W_GR_dT);
        ip_cv.dfW_2b_dT =
            drho_W_FR_dT * pCap * (ip_cv.biot_data() - ip_data.phi) *
                ip_cv.beta_p_SR()
#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
            - rho_W_FR * pCap * ip_data.dphi_dT * ip_cv.beta_p_SR()
#endif
            ;

        if (dt == 0.)
        {
            ip_cv.dfW_3a_dp_GR = 0.;
            ip_cv.dfW_3a_dp_cap = 0.;
            ip_cv.dfW_3a_dT = 0.;
        }
        else
        {
            double const rho_W_GR_dot =
                (ip_data.rhoWGR - ip_data.rhoWGR_prev) / dt;
            double const rho_W_LR_dot =
                (ip_data.rhoWLR - ip_data.rhoWLR_prev) / dt;

            ip_cv.dfW_3a_dp_GR =
                /*(ds_G_dp_GR = 0) * rho_W_GR_dot +*/ s_G * c.drho_W_GR_dp_GR /
                    dt +
                /*(ds_L_dp_GR = 0) * rho_W_LR_dot +*/ s_L * c.drho_W_LR_dp_GR /
                    dt;
            ip_cv.dfW_3a_dp_cap = ds_G_dp_cap * rho_W_GR_dot +
                                  s_G * c.drho_W_GR_dp_cap / dt +
                                  ip_cv.dS_L_dp_cap() * rho_W_LR_dot -
                                  s_L * c.drho_W_LR_dp_LR / dt;
            ip_cv.dfW_3a_dT =
                s_G * c.drho_W_GR_dT / dt + s_L * c.drho_W_LR_dT / dt;
        }

        ip_cv.dfW_4_LWpG_a_dp_GR = c.drho_W_GR_dp_GR * k_over_mu_G
                                   // + rhoWGR * (dk_over_mu_G_dp_GR = 0)
                                   + c.drho_W_LR_dp_GR * k_over_mu_L
            // + rhoWLR * (dk_over_mu_L_dp_GR = 0)
            ;
        ip_cv.dfW_4_LWpG_a_dp_cap = c.drho_W_GR_dp_cap * k_over_mu_G +
                                    ip_data.rhoWGR * ip_cv.dk_over_mu_G_dp_cap +
                                    -c.drho_W_LR_dp_LR * k_over_mu_L +
                                    ip_data.rhoWLR * ip_cv.dk_over_mu_L_dp_cap;

        ip_cv.dfW_4_LWpG_a_dT =
            c.drho_W_GR_dT * k_over_mu_G
            //+ rhoWGR * (dk_over_mu_G_dT != 0 TODO for mu_G(T))
            + c.drho_W_LR_dT * k_over_mu_L
            //+ rhoWLR * (dk_over_mu_L_dT != 0 TODO for mu_G(T))
            ;

        // TODO (naumov) for dxmW*/d* != 0
        ip_cv.dfW_4_LWpG_d_dp_GR =
            Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Zero();
        ip_cv.dfW_4_LWpG_d_dp_cap =
            Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Zero();
        ip_cv.dfW_4_LWpG_d_dT =
            Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Zero();

        ip_cv.dfW_4_LWpC_a_dp_GR = c.drho_W_LR_dp_GR * k_over_mu_L
            //+ rhoWLR * (dk_over_mu_L_dp_GR = 0)
            ;
        ip_cv.dfW_4_LWpC_a_dp_cap = -c.drho_W_LR_dp_LR * k_over_mu_L +
                                    ip_data.rhoWLR * ip_cv.dk_over_mu_L_dp_cap;
        ip_cv.dfW_4_LWpC_a_dT = c.drho_W_LR_dT * k_over_mu_L
            //+ rhoWLR * (dk_over_mu_L_dT != 0 TODO for mu_L(T))
            ;

        // TODO (naumov) for dxmW*/d* != 0
        ip_cv.dfW_4_LWpC_d_dp_GR =
            Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Zero();
        ip_cv.dfW_4_LWpC_d_dp_cap =
            Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Zero();
        ip_cv.dfW_4_LWpC_d_dT =
            Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Zero();

        ip_cv.dfC_4_LCpC_a_dp_GR = c.drho_C_LR_dp_GR * k_over_mu_L
            //+ rhoCLR * (dk_over_mu_L_dp_GR = 0)
            ;
        ip_cv.dfC_4_LCpC_a_dp_cap = -c.drho_C_LR_dp_LR * k_over_mu_L +
                                    ip_data.rhoCLR * ip_cv.dk_over_mu_L_dp_cap;
        ip_cv.dfC_4_LCpC_a_dT = c.drho_W_LR_dT * k_over_mu_L
            //+ rhoWLR * (dk_over_mu_L_dT != 0 TODO for mu_L(T))
            ;

        // TODO (naumov) for dxmW*/d* != 0
        ip_cv.dfC_4_LCpC_d_dp_GR =
            Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Zero();
        ip_cv.dfC_4_LCpC_d_dp_cap =
            Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Zero();
        ip_cv.dfC_4_LCpC_d_dT =
            Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Zero();
    }

    return {ip_constitutive_data, ip_constitutive_variables};
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::size_t TH2MLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::setIPDataInitialConditions(std::string_view name,
                                                 double const* values,
                                                 int const integration_order)
{
    if (integration_order !=
        static_cast<int>(this->integration_method_.getIntegrationOrder()))
    {
        OGS_FATAL(
            "Setting integration point initial conditions; The integration "
            "order of the local assembler for element {:d} is different "
            "from the integration order in the initial condition.",
            this->element_.getID());
    }

    if (name == "sigma" && this->process_data_.initial_stress.value)
    {
        OGS_FATAL(
            "Setting initial conditions for stress from integration "
            "point data and from a parameter '{:s}' is not possible "
            "simultaneously.",
            this->process_data_.initial_stress.value->name);
    }

    if (name.starts_with("material_state_variable_"))
    {
        name.remove_prefix(24);
        DBUG("Setting material state variable '{:s}'", name);

        auto const& internal_variables =
            this->solid_material_.getInternalVariables();
        if (auto const iv = std::find_if(
                begin(internal_variables), end(internal_variables),
                [&name](auto const& iv) { return iv.name == name; });
            iv != end(internal_variables))
        {
            DBUG("Setting material state variable '{:s}'", name);
            return ProcessLib::setIntegrationPointDataMaterialStateVariables(
                values, this->material_states_,
                &ConstitutiveRelations::MaterialStateData<
                    DisplacementDim>::material_state_variables,
                iv->reference);
        }

        WARN(
            "Could not find variable {:s} in solid material model's "
            "internal variables.",
            name);
        return 0;
    }

    // TODO this logic could be pulled out of the local assembler into the
    // process. That might lead to a slightly better performance due to less
    // string comparisons.
    return ProcessLib::Reflection::reflectSetIPData<DisplacementDim>(
        name, values, this->current_states_);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void TH2MLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                        DisplacementDim>::
    setInitialConditionsConcrete(Eigen::VectorXd const local_x,
                                 double const t,
                                 int const /*process_id*/)
{
    [[maybe_unused]] auto const matrix_size =
        gas_pressure_size + capillary_pressure_size + temperature_size +
        displacement_size;

    assert(local_x.size() == matrix_size);

    auto const capillary_pressure =
        local_x.template segment<capillary_pressure_size>(
            capillary_pressure_index);

    auto const p_GR =
        local_x.template segment<gas_pressure_size>(gas_pressure_index);

    auto const temperature =
        local_x.template segment<temperature_size>(temperature_index);

    auto const displacement =
        local_x.template segment<displacement_size>(displacement_index);

    constexpr double dt = std::numeric_limits<double>::quiet_NaN();
    auto const& medium =
        *this->process_data_.media_map.getMedium(this->element_.getID());
    auto const& solid_phase = medium.phase("Solid");

    ConstitutiveRelations::ConstitutiveModels<DisplacementDim> const models{
        this->solid_material_};

    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        MPL::VariableArray vars;

        auto& ip_data = _ip_data[ip];
        auto& ip_out = this->output_data_[ip];
        auto& prev_state = this->prev_states_[ip];
        auto const& Np = ip_data.N_p;
        auto const& NT = Np;
        auto const& Nu = ip_data.N_u;
        auto const& gradNu = ip_data.dNdx_u;
        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                this->element_, Nu);
        ParameterLib::SpatialPosition const pos{
            std::nullopt, this->element_.getID(), ip,
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    this->element_, ip_data.N_u))};

        double const pCap = Np.dot(capillary_pressure);
        vars.capillary_pressure = pCap;

        double const T = NT.dot(temperature);
        TemperatureData const T_data{T, T};  // T_prev = T in initialization.
        vars.temperature = T;

        auto const Bu =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                gradNu, Nu, x_coord, this->is_axially_symmetric_);

        auto& eps = ip_out.eps_data.eps;
        eps.noalias() = Bu * displacement;

        // Set volumetric strain rate for the general case without swelling.
        vars.volumetric_strain = Invariants::trace(eps);

        double const S_L =
            medium.property(MPL::PropertyType::saturation)
                .template value<double>(
                    vars, pos, t, std::numeric_limits<double>::quiet_NaN());
        this->prev_states_[ip].S_L_data->S_L = S_L;

        // TODO (naumov) Double computation of C_el might be avoided if
        // updateConstitutiveVariables is called before. But it might interfere
        // with eps_m initialization.
        ConstitutiveRelations::ElasticTangentStiffnessData<DisplacementDim>
            C_el_data;
        models.elastic_tangent_stiffness_model.eval({pos, t, dt}, T_data,
                                                    C_el_data);
        auto const& C_el = C_el_data.stiffness_tensor;

        // Set eps_m_prev from potentially non-zero eps and sigma_sw from
        // restart.
        auto const& sigma_sw = this->current_states_[ip].swelling_data.sigma_sw;
        prev_state.mechanical_strain_data->eps_m.noalias() =
            solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate)
                ? eps + C_el.inverse() * sigma_sw
                : eps;

        if (this->process_data_.initial_stress.isTotalStress())
        {
            auto const alpha_b =
                medium.property(MPL::PropertyType::biot_coefficient)
                    .template value<double>(vars, pos, t, 0.0 /*dt*/);

            vars.liquid_saturation = S_L;
            double const bishop =
                medium.property(MPL::PropertyType::bishops_effective_stress)
                    .template value<double>(vars, pos, t, 0.0 /*dt*/);

            this->current_states_[ip].eff_stress_data.sigma.noalias() +=
                alpha_b * Np.dot(p_GR - bishop * capillary_pressure) *
                Invariants::identity2;
            this->prev_states_[ip].eff_stress_data =
                this->current_states_[ip].eff_stress_data;
        }
    }

    // local_x_prev equal to local_x s.t. the local_x_dot is zero.
    updateConstitutiveVariables(local_x, local_x, t, 0, models);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto& ip_data = _ip_data[ip];
        ip_data.pushBackState();

        this->material_states_[ip].pushBackState();
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void TH2MLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::assemble(double const t, double const dt,
                               std::vector<double> const& local_x,
                               std::vector<double> const& local_x_prev,
                               std::vector<double>& local_M_data,
                               std::vector<double>& local_K_data,
                               std::vector<double>& local_rhs_data)
{
    auto const matrix_size = gas_pressure_size + capillary_pressure_size +
                             temperature_size + displacement_size;
    assert(local_x.size() == matrix_size);

    auto const capillary_pressure =
        Eigen::Map<VectorType<capillary_pressure_size> const>(
            local_x.data() + capillary_pressure_index, capillary_pressure_size);

    auto const capillary_pressure_prev =
        Eigen::Map<VectorType<capillary_pressure_size> const>(
            local_x_prev.data() + capillary_pressure_index,
            capillary_pressure_size);

    // pointer to local_M_data vector
    auto local_M =
        MathLib::createZeroedMatrix<MatrixType<matrix_size, matrix_size>>(
            local_M_data, matrix_size, matrix_size);

    // pointer to local_K_data vector
    auto local_K =
        MathLib::createZeroedMatrix<MatrixType<matrix_size, matrix_size>>(
            local_K_data, matrix_size, matrix_size);

    // pointer to local_rhs_data vector
    auto local_f = MathLib::createZeroedVector<VectorType<matrix_size>>(
        local_rhs_data, matrix_size);

    // component-formulation
    // W - liquid phase main component
    // C - gas phase main component
    // pointer-matrices to the mass matrix - C component equation
    auto MCpG = local_M.template block<C_size, gas_pressure_size>(
        C_index, gas_pressure_index);
    auto MCpC = local_M.template block<C_size, capillary_pressure_size>(
        C_index, capillary_pressure_index);
    auto MCT = local_M.template block<C_size, temperature_size>(
        C_index, temperature_index);
    auto MCu = local_M.template block<C_size, displacement_size>(
        C_index, displacement_index);

    // pointer-matrices to the stiffness matrix - C component equation
    auto LCpG = local_K.template block<C_size, gas_pressure_size>(
        C_index, gas_pressure_index);
    auto LCpC = local_K.template block<C_size, capillary_pressure_size>(
        C_index, capillary_pressure_index);
    auto LCT = local_K.template block<C_size, temperature_size>(
        C_index, temperature_index);

    // pointer-matrices to the mass matrix - W component equation
    auto MWpG = local_M.template block<W_size, gas_pressure_size>(
        W_index, gas_pressure_index);
    auto MWpC = local_M.template block<W_size, capillary_pressure_size>(
        W_index, capillary_pressure_index);
    auto MWT = local_M.template block<W_size, temperature_size>(
        W_index, temperature_index);
    auto MWu = local_M.template block<W_size, displacement_size>(
        W_index, displacement_index);

    // pointer-matrices to the stiffness matrix - W component equation
    auto LWpG = local_K.template block<W_size, gas_pressure_size>(
        W_index, gas_pressure_index);
    auto LWpC = local_K.template block<W_size, capillary_pressure_size>(
        W_index, capillary_pressure_index);
    auto LWT = local_K.template block<W_size, temperature_size>(
        W_index, temperature_index);

    // pointer-matrices to the mass matrix - temperature equation
    auto MTu = local_M.template block<temperature_size, displacement_size>(
        temperature_index, displacement_index);

    // pointer-matrices to the stiffness matrix - temperature equation
    auto KTT = local_K.template block<temperature_size, temperature_size>(
        temperature_index, temperature_index);

    // pointer-matrices to the stiffness matrix - displacement equation
    auto KUpG = local_K.template block<displacement_size, gas_pressure_size>(
        displacement_index, gas_pressure_index);
    auto KUpC =
        local_K.template block<displacement_size, capillary_pressure_size>(
            displacement_index, capillary_pressure_index);

    // pointer-vectors to the right hand side terms - C-component equation
    auto fC = local_f.template segment<C_size>(C_index);
    // pointer-vectors to the right hand side terms - W-component equation
    auto fW = local_f.template segment<W_size>(W_index);
    // pointer-vectors to the right hand side terms - temperature equation
    auto fT = local_f.template segment<temperature_size>(temperature_index);
    // pointer-vectors to the right hand side terms - displacement equation
    auto fU = local_f.template segment<displacement_size>(displacement_index);

    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();

    ConstitutiveRelations::ConstitutiveModels<DisplacementDim> const models{
        this->solid_material_};

    auto const [ip_constitutive_data, ip_constitutive_variables] =
        updateConstitutiveVariables(
            Eigen::Map<Eigen::VectorXd const>(local_x.data(), local_x.size()),
            Eigen::Map<Eigen::VectorXd const>(local_x_prev.data(),
                                              local_x_prev.size()),
            t, dt, models);

    for (unsigned int_point = 0; int_point < n_integration_points; int_point++)
    {
        auto& ip = _ip_data[int_point];
        auto& ip_cv = ip_constitutive_variables[int_point];
        auto& ip_out = this->output_data_[int_point];
        auto& current_state = this->current_states_[int_point];
        auto const& prev_state = this->prev_states_[int_point];

        auto const& Np = ip.N_p;
        auto const& NT = Np;
        auto const& Nu = ip.N_u;
        ParameterLib::SpatialPosition const pos{
            std::nullopt, this->element_.getID(), int_point,
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    this->element_, Nu))};

        auto const& NpT = Np.transpose().eval();
        auto const& NTT = NT.transpose().eval();

        auto const& gradNp = ip.dNdx_p;
        auto const& gradNT = gradNp;
        auto const& gradNu = ip.dNdx_u;

        auto const& gradNpT = gradNp.transpose().eval();
        auto const& gradNTT = gradNT.transpose().eval();

        auto const& w = ip.integration_weight;

        auto const& m = Invariants::identity2;

        auto const mT = m.transpose().eval();

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                this->element_, Nu);

        auto const Bu =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                gradNu, Nu, x_coord, this->is_axially_symmetric_);

        auto const BuT = Bu.transpose().eval();

        double const pCap = Np.dot(capillary_pressure);
        double const pCap_prev = Np.dot(capillary_pressure_prev);

        double const beta_T_SR = ip_cv.s_therm_exp_data.beta_T_SR;

        auto const I =
            Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Identity();

        const double sD_G = ip.diffusion_coefficient_vapour;
        const double sD_L = ip.diffusion_coefficient_solute;

        GlobalDimMatrixType const D_C_G = sD_G * I;
        GlobalDimMatrixType const D_W_G = sD_G * I;
        GlobalDimMatrixType const D_C_L = sD_L * I;
        GlobalDimMatrixType const D_W_L = sD_L * I;

        auto const s_L = current_state.S_L_data.S_L;
        auto const s_G = 1. - s_L;
        auto const s_L_dot = (s_L - prev_state.S_L_data->S_L) / dt;

        auto& alpha_B = ip_cv.biot_data();
        auto& beta_p_SR = ip_cv.beta_p_SR();

        auto const& b = this->process_data_.specific_body_force;

        // porosity
        auto& phi = ip.phi;

        // volume fraction
        auto const phi_G = s_G * phi;
        auto const phi_L = s_L * phi;
        auto const phi_S = 1. - phi;

        // solid phase density
        auto& rho_SR = ip.rhoSR;
        // effective density
        auto const rho = phi_G * ip.rhoGR + phi_L * ip.rhoLR + phi_S * rho_SR;

        // abbreviations
        const double rho_C_FR = s_G * ip.rhoCGR + s_L * ip.rhoCLR;
        const double rho_W_FR = s_G * ip.rhoWGR + s_L * ip.rhoWLR;

        // phase specific enthalpies
        auto& h_G = ip.h_G;
        auto& h_L = ip.h_L;

        auto const rho_C_GR_dot = (ip.rhoCGR - ip.rhoCGR_prev) / dt;
        auto const rho_C_LR_dot = (ip.rhoCLR - ip.rhoCLR_prev) / dt;
        auto const rho_W_GR_dot = (ip.rhoWGR - ip.rhoWGR_prev) / dt;
        auto const rho_W_LR_dot = (ip.rhoWLR - ip.rhoWLR_prev) / dt;

        auto const rho_h_eff = ip.rho_G_h_G + ip.rho_L_h_L + ip.rho_S_h_S;

        auto const rho_u_eff_dot = (ip.rho_u_eff - ip.rho_u_eff_prev) / dt;

        GlobalDimMatrixType const k_over_mu_G =
            ip_out.permeability_data.Ki * ip_out.permeability_data.k_rel_G /
            ip.muGR;
        GlobalDimMatrixType const k_over_mu_L =
            ip_out.permeability_data.Ki * ip_out.permeability_data.k_rel_L /
            ip.muLR;

        // ---------------------------------------------------------------------
        // C-component equation
        // ---------------------------------------------------------------------

        MCpG.noalias() += NpT * rho_C_FR * (alpha_B - phi) * beta_p_SR * Np * w;
        MCpC.noalias() -=
            NpT * rho_C_FR * (alpha_B - phi) * beta_p_SR * s_L * Np * w;

        if (this->process_data_.apply_mass_lumping)
        {
            if (pCap - pCap_prev != 0.)  // avoid division by Zero
            {
                MCpC.noalias() +=
                    NpT *
                    (phi * (ip.rhoCLR - ip.rhoCGR) -
                     rho_C_FR * pCap * (alpha_B - phi) * beta_p_SR) *
                    s_L_dot * dt / (pCap - pCap_prev) * Np * w;
            }
        }

        MCT.noalias() -= NpT * rho_C_FR * (alpha_B - phi) * beta_T_SR * Np * w;
        MCu.noalias() += NpT * rho_C_FR * alpha_B * mT * Bu * w;

        using DisplacementDimMatrix =
            Eigen::Matrix<double, DisplacementDim, DisplacementDim>;

        DisplacementDimMatrix const advection_C_G = ip.rhoCGR * k_over_mu_G;
        DisplacementDimMatrix const advection_C_L = ip.rhoCLR * k_over_mu_L;

        DisplacementDimMatrix const diffusion_CGpGR =
            -phi_G * ip.rhoGR * D_C_G * ip.dxmWG_dpGR;
        DisplacementDimMatrix const diffusion_CLpGR =
            -phi_L * ip.rhoLR * D_C_L * ip.dxmWL_dpGR;

        DisplacementDimMatrix const diffusion_CGpCap =
            -phi_G * ip.rhoGR * D_C_G * ip.dxmWG_dpCap;
        DisplacementDimMatrix const diffusion_CLpCap =
            -phi_L * ip.rhoLR * D_C_L * ip.dxmWL_dpCap;

        DisplacementDimMatrix const diffusion_CGT =
            -phi_G * ip.rhoGR * D_C_G * ip.dxmWG_dT;
        DisplacementDimMatrix const diffusion_CLT =
            -phi_L * ip.rhoLR * D_C_L * ip.dxmWL_dT;

        DisplacementDimMatrix const advection_C = advection_C_G + advection_C_L;
        DisplacementDimMatrix const diffusion_C_pGR =
            diffusion_CGpGR + diffusion_CLpGR;
        DisplacementDimMatrix const diffusion_C_pCap =
            diffusion_CGpCap + diffusion_CLpCap;

        DisplacementDimMatrix const diffusion_C_T =
            diffusion_CGT + diffusion_CLT;

        LCpG.noalias() +=
            gradNpT * (advection_C + diffusion_C_pGR) * gradNp * w;

        LCpC.noalias() +=
            gradNpT * (diffusion_C_pCap - advection_C_L) * gradNp * w;

        LCT.noalias() += gradNpT * (diffusion_C_T)*gradNp * w;

        fC.noalias() += gradNpT *
                        (advection_C_G * ip.rhoGR + advection_C_L * ip.rhoLR) *
                        b * w;

        if (!this->process_data_.apply_mass_lumping)
        {
            fC.noalias() -= NpT *
                            (phi * (ip.rhoCLR - ip.rhoCGR) -
                             rho_C_FR * pCap * (alpha_B - phi) * beta_p_SR) *
                            s_L_dot * w;
        }
        // fC_III
        fC.noalias() -=
            NpT * phi * (s_G * rho_C_GR_dot + s_L * rho_C_LR_dot) * w;

        // ---------------------------------------------------------------------
        // W-component equation
        // ---------------------------------------------------------------------

        MWpG.noalias() += NpT * rho_W_FR * (alpha_B - phi) * beta_p_SR * Np * w;
        MWpC.noalias() -=
            NpT * rho_W_FR * (alpha_B - phi) * beta_p_SR * s_L * Np * w;

        if (this->process_data_.apply_mass_lumping)
        {
            if (pCap - pCap_prev != 0.)  // avoid division by Zero
            {
                MWpC.noalias() +=
                    NpT *
                    (phi * (ip.rhoWLR - ip.rhoWGR) -
                     rho_W_FR * pCap * (alpha_B - phi) * beta_p_SR) *
                    s_L_dot * dt / (pCap - pCap_prev) * Np * w;
            }
        }

        MWT.noalias() -= NpT * rho_W_FR * (alpha_B - phi) * beta_T_SR * Np * w;

        MWu.noalias() += NpT * rho_W_FR * alpha_B * mT * Bu * w;

        DisplacementDimMatrix const advection_W_G = ip.rhoWGR * k_over_mu_G;
        DisplacementDimMatrix const advection_W_L = ip.rhoWLR * k_over_mu_L;

        DisplacementDimMatrix const diffusion_WGpGR =
            phi_G * ip.rhoGR * D_W_G * ip.dxmWG_dpGR;
        DisplacementDimMatrix const diffusion_WLpGR =
            phi_L * ip.rhoLR * D_W_L * ip.dxmWL_dpGR;

        DisplacementDimMatrix const diffusion_WGpCap =
            phi_G * ip.rhoGR * D_W_G * ip.dxmWG_dpCap;
        DisplacementDimMatrix const diffusion_WLpCap =
            phi_L * ip.rhoLR * D_W_L * ip.dxmWL_dpCap;

        DisplacementDimMatrix const diffusion_WGT =
            phi_G * ip.rhoGR * D_W_G * ip.dxmWG_dT;
        DisplacementDimMatrix const diffusion_WLT =
            phi_L * ip.rhoLR * D_W_L * ip.dxmWL_dT;

        DisplacementDimMatrix const advection_W = advection_W_G + advection_W_L;
        DisplacementDimMatrix const diffusion_W_pGR =
            diffusion_WGpGR + diffusion_WLpGR;
        DisplacementDimMatrix const diffusion_W_pCap =
            diffusion_WGpCap + diffusion_WLpCap;

        DisplacementDimMatrix const diffusion_W_T =
            diffusion_WGT + diffusion_WLT;

        LWpG.noalias() +=
            gradNpT * (advection_W + diffusion_W_pGR) * gradNp * w;

        LWpC.noalias() +=
            gradNpT * (diffusion_W_pCap - advection_W_L) * gradNp * w;

        LWT.noalias() += gradNpT * (diffusion_W_T)*gradNp * w;

        fW.noalias() += gradNpT *
                        (advection_W_G * ip.rhoGR + advection_W_L * ip.rhoLR) *
                        b * w;

        if (!this->process_data_.apply_mass_lumping)
        {
            fW.noalias() -= NpT *
                            (phi * (ip.rhoWLR - ip.rhoWGR) -
                             rho_W_FR * pCap * (alpha_B - phi) * beta_p_SR) *
                            s_L_dot * w;
        }

        fW.noalias() -=
            NpT * phi * (s_G * rho_W_GR_dot + s_L * rho_W_LR_dot) * w;

        // ---------------------------------------------------------------------
        //  - temperature equation
        // ---------------------------------------------------------------------

        MTu.noalias() += NTT * rho_h_eff * mT * Bu * w;

        KTT.noalias() += gradNTT * ip.lambda * gradNT * w;

        fT.noalias() -= NTT * rho_u_eff_dot * w;

        fT.noalias() +=
            gradNTT * (ip.rhoGR * h_G * ip.w_GS + ip.rhoLR * h_L * ip.w_LS) * w;

        fT.noalias() +=
            gradNTT *
            (ip.rhoCGR * ip.h_CG * ip.d_CG + ip.rhoWGR * ip.h_WG * ip.d_WG) * w;

        fT.noalias() +=
            NTT *
            (ip.rhoGR * ip.w_GS.transpose() + ip.rhoLR * ip.w_LS.transpose()) *
            b * w;

        // ---------------------------------------------------------------------
        //  - displacement equation
        // ---------------------------------------------------------------------

        KUpG.noalias() -= (BuT * alpha_B * m * Np) * w;

        KUpC.noalias() += (BuT * alpha_B * ip_cv.chi_S_L.chi_S_L * m * Np) * w;

        fU.noalias() -= (BuT * current_state.eff_stress_data.sigma -
                         N_u_op(Nu).transpose() * rho * b) *
                        w;

        if (this->process_data_.apply_mass_lumping)
        {
            MCpG = MCpG.colwise().sum().eval().asDiagonal();
            MCpC = MCpC.colwise().sum().eval().asDiagonal();
            MWpG = MWpG.colwise().sum().eval().asDiagonal();
            MWpC = MWpC.colwise().sum().eval().asDiagonal();
        }
    }  // int_point-loop
}

// Assembles the local Jacobian matrix. So far, the linearisation of HT part is
// not considered as that in HT process.
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void TH2MLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                        DisplacementDim>::
    assembleWithJacobian(double const t, double const dt,
                         std::vector<double> const& local_x,
                         std::vector<double> const& local_x_prev,
                         std::vector<double>& /*local_M_data*/,
                         std::vector<double>& /*local_K_data*/,
                         std::vector<double>& local_rhs_data,
                         std::vector<double>& local_Jac_data)
{
    auto const matrix_size = gas_pressure_size + capillary_pressure_size +
                             temperature_size + displacement_size;
    assert(local_x.size() == matrix_size);

    auto const temperature = Eigen::Map<VectorType<temperature_size> const>(
        local_x.data() + temperature_index, temperature_size);

    auto const gas_pressure = Eigen::Map<VectorType<gas_pressure_size> const>(
        local_x.data() + gas_pressure_index, gas_pressure_size);

    auto const capillary_pressure =
        Eigen::Map<VectorType<capillary_pressure_size> const>(
            local_x.data() + capillary_pressure_index, capillary_pressure_size);

    auto const displacement = Eigen::Map<VectorType<displacement_size> const>(
        local_x.data() + displacement_index, displacement_size);

    auto const gas_pressure_prev =
        Eigen::Map<VectorType<gas_pressure_size> const>(
            local_x_prev.data() + gas_pressure_index, gas_pressure_size);

    auto const capillary_pressure_prev =
        Eigen::Map<VectorType<capillary_pressure_size> const>(
            local_x_prev.data() + capillary_pressure_index,
            capillary_pressure_size);

    auto const temperature_prev =
        Eigen::Map<VectorType<temperature_size> const>(
            local_x_prev.data() + temperature_index, temperature_size);

    auto const displacement_prev =
        Eigen::Map<VectorType<displacement_size> const>(
            local_x_prev.data() + displacement_index, displacement_size);

    auto local_Jac =
        MathLib::createZeroedMatrix<MatrixType<matrix_size, matrix_size>>(
            local_Jac_data, matrix_size, matrix_size);

    auto local_f = MathLib::createZeroedVector<VectorType<matrix_size>>(
        local_rhs_data, matrix_size);

    // component-formulation
    // W - liquid phase main component
    // C - gas phase main component

    // C component equation matrices
    MatrixType<C_size, gas_pressure_size> MCpG =
        MatrixType<C_size, gas_pressure_size>::Zero(C_size, gas_pressure_size);
    MatrixType<C_size, capillary_pressure_size> MCpC =
        MatrixType<C_size, capillary_pressure_size>::Zero(
            C_size, capillary_pressure_size);
    MatrixType<C_size, temperature_size> MCT =
        MatrixType<C_size, temperature_size>::Zero(C_size, temperature_size);
    MatrixType<C_size, displacement_size> MCu =
        MatrixType<C_size, displacement_size>::Zero(C_size, displacement_size);

    MatrixType<C_size, gas_pressure_size> LCpG =
        MatrixType<C_size, gas_pressure_size>::Zero(C_size, gas_pressure_size);
    MatrixType<C_size, capillary_pressure_size> LCpC =
        MatrixType<C_size, capillary_pressure_size>::Zero(
            C_size, capillary_pressure_size);
    MatrixType<C_size, temperature_size> LCT =
        MatrixType<C_size, temperature_size>::Zero(C_size, temperature_size);

    // mass matrix - W component equation
    MatrixType<W_size, gas_pressure_size> MWpG =
        MatrixType<W_size, gas_pressure_size>::Zero(W_size, gas_pressure_size);
    MatrixType<W_size, capillary_pressure_size> MWpC =
        MatrixType<W_size, capillary_pressure_size>::Zero(
            W_size, capillary_pressure_size);
    MatrixType<W_size, temperature_size> MWT =
        MatrixType<W_size, temperature_size>::Zero(W_size, temperature_size);
    MatrixType<W_size, displacement_size> MWu =
        MatrixType<W_size, displacement_size>::Zero(W_size, displacement_size);

    // stiffness matrix - W component equation
    MatrixType<W_size, gas_pressure_size> LWpG =
        MatrixType<W_size, gas_pressure_size>::Zero(W_size, gas_pressure_size);
    MatrixType<W_size, capillary_pressure_size> LWpC =
        MatrixType<W_size, capillary_pressure_size>::Zero(
            W_size, capillary_pressure_size);
    MatrixType<W_size, temperature_size> LWT =
        MatrixType<W_size, temperature_size>::Zero(W_size, temperature_size);

    // mass matrix - temperature equation
    MatrixType<temperature_size, displacement_size> MTu =
        MatrixType<temperature_size, displacement_size>::Zero(
            temperature_size, displacement_size);

    // stiffness matrix - temperature equation
    MatrixType<temperature_size, temperature_size> KTT =
        MatrixType<temperature_size, temperature_size>::Zero(temperature_size,
                                                             temperature_size);

    // stiffness matrices - displacement equation coupling into pressures
    MatrixType<displacement_size, gas_pressure_size> KUpG =
        MatrixType<displacement_size, gas_pressure_size>::Zero(
            displacement_size, gas_pressure_size);
    MatrixType<displacement_size, capillary_pressure_size> KUpC =
        MatrixType<displacement_size, capillary_pressure_size>::Zero(
            displacement_size, capillary_pressure_size);

    // pointer-vectors to the right hand side terms - C-component equation
    auto fC = local_f.template segment<C_size>(C_index);
    // pointer-vectors to the right hand side terms - W-component equation
    auto fW = local_f.template segment<W_size>(W_index);
    // pointer-vectors to the right hand side terms - temperature equation
    auto fT = local_f.template segment<temperature_size>(temperature_index);
    // pointer-vectors to the right hand side terms - displacement equation
    auto fU = local_f.template segment<displacement_size>(displacement_index);

    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();

    ConstitutiveRelations::ConstitutiveModels<DisplacementDim> const models{
        this->solid_material_};

    auto const [ip_constitutive_data, ip_constitutive_variables] =
        updateConstitutiveVariables(
            Eigen::Map<Eigen::VectorXd const>(local_x.data(), local_x.size()),
            Eigen::Map<Eigen::VectorXd const>(local_x_prev.data(),
                                              local_x_prev.size()),
            t, dt, models);

    for (unsigned int_point = 0; int_point < n_integration_points; int_point++)
    {
        auto& ip = _ip_data[int_point];
        auto& ip_cd = ip_constitutive_data[int_point];
        auto& ip_cv = ip_constitutive_variables[int_point];
        auto& ip_out = this->output_data_[int_point];
        auto& current_state = this->current_states_[int_point];
        auto& prev_state = this->prev_states_[int_point];

        auto const& Np = ip.N_p;
        auto const& NT = Np;
        auto const& Nu = ip.N_u;
        ParameterLib::SpatialPosition const pos{
            std::nullopt, this->element_.getID(), int_point,
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    this->element_, Nu))};

        auto const& NpT = Np.transpose().eval();
        auto const& NTT = NT.transpose().eval();

        auto const& gradNp = ip.dNdx_p;
        auto const& gradNT = gradNp;
        auto const& gradNu = ip.dNdx_u;

        auto const& gradNpT = gradNp.transpose().eval();
        auto const& gradNTT = gradNT.transpose().eval();

        auto const& w = ip.integration_weight;

        auto const& m = Invariants::identity2;
        auto const mT = m.transpose().eval();

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                this->element_, Nu);

        auto const Bu =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                gradNu, Nu, x_coord, this->is_axially_symmetric_);

        auto const BuT = Bu.transpose().eval();

        double const div_u_dot =
            Invariants::trace(Bu * (displacement - displacement_prev) / dt);

        double const pGR = Np.dot(gas_pressure);
        double const pCap = Np.dot(capillary_pressure);
        double const T = NT.dot(temperature);

        GlobalDimVectorType const gradpGR = gradNp * gas_pressure;
        GlobalDimVectorType const gradpCap = gradNp * capillary_pressure;
        GlobalDimVectorType const gradT = gradNT * temperature;

        double const pGR_prev = Np.dot(gas_pressure_prev);
        double const pCap_prev = Np.dot(capillary_pressure_prev);
        double const T_prev = NT.dot(temperature_prev);
        double const beta_T_SR = ip_cv.s_therm_exp_data.beta_T_SR;

        auto const I =
            Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Identity();

        const double sD_G = ip.diffusion_coefficient_vapour;
        const double sD_L = ip.diffusion_coefficient_solute;

        GlobalDimMatrixType const D_C_G = sD_G * I;
        GlobalDimMatrixType const D_W_G = sD_G * I;
        GlobalDimMatrixType const D_C_L = sD_L * I;
        GlobalDimMatrixType const D_W_L = sD_L * I;

        auto& s_L = current_state.S_L_data.S_L;
        auto const s_G = 1. - s_L;
        auto const s_L_dot = (s_L - prev_state.S_L_data->S_L) / dt;

        auto const alpha_B = ip_cv.biot_data();
        auto const beta_p_SR = ip_cv.beta_p_SR();

        auto const& b = this->process_data_.specific_body_force;

        // porosity
        auto& phi = ip.phi;

        // volume fraction
        auto const phi_G = s_G * phi;
        auto const phi_L = s_L * phi;
        auto const phi_S = 1. - phi;

        // solid phase density
        auto& rho_SR = ip.rhoSR;
        // effective density
        auto const rho = phi_G * ip.rhoGR + phi_L * ip.rhoLR + phi_S * rho_SR;

        // abbreviations
        const double rho_C_FR = s_G * ip.rhoCGR + s_L * ip.rhoCLR;
        const double rho_W_FR = s_G * ip.rhoWGR + s_L * ip.rhoWLR;

        // phase specific enthalpies
        auto& h_G = ip.h_G;
        auto& h_L = ip.h_L;

        auto const rho_C_GR_dot = (ip.rhoCGR - ip.rhoCGR_prev) / dt;
        auto const rho_C_LR_dot = (ip.rhoCLR - ip.rhoCLR_prev) / dt;
        auto const rho_W_GR_dot = (ip.rhoWGR - ip.rhoWGR_prev) / dt;
        auto const rho_W_LR_dot = (ip.rhoWLR - ip.rhoWLR_prev) / dt;

        auto const rho_h_eff = ip.rho_G_h_G + ip.rho_L_h_L + ip.rho_S_h_S;

        auto const rho_u_eff_dot = (ip.rho_u_eff - ip.rho_u_eff_prev) / dt;

        GlobalDimMatrixType const k_over_mu_G =
            ip_out.permeability_data.Ki * ip_out.permeability_data.k_rel_G /
            ip.muGR;
        GlobalDimMatrixType const k_over_mu_L =
            ip_out.permeability_data.Ki * ip_out.permeability_data.k_rel_L /
            ip.muLR;

        // ---------------------------------------------------------------------
        // C-component equation
        // ---------------------------------------------------------------------

        MCpG.noalias() += NpT * rho_C_FR * (alpha_B - phi) * beta_p_SR * Np * w;
        MCpC.noalias() -=
            NpT * rho_C_FR * (alpha_B - phi) * beta_p_SR * s_L * Np * w;

        if (this->process_data_.apply_mass_lumping)
        {
            if (pCap - pCap_prev != 0.)  // avoid division by Zero
            {
                MCpC.noalias() +=
                    NpT *
                    (phi * (ip.rhoCLR - ip.rhoCGR) -
                     rho_C_FR * pCap * (alpha_B - phi) * beta_p_SR) *
                    s_L_dot * dt / (pCap - pCap_prev) * Np * w;
            }
        }

        MCT.noalias() -= NpT * rho_C_FR * (alpha_B - phi) * beta_T_SR * Np * w;
        // d (fC_4_MCT * T_dot)/d T
        local_Jac
            .template block<C_size, temperature_size>(C_index,
                                                      temperature_index)
            .noalias() += NpT * ip_cv.dfC_4_MCT_dT * (T - T_prev) / dt * NT * w;

        MCu.noalias() += NpT * rho_C_FR * alpha_B * mT * Bu * w;
        // d (fC_4_MCu * u_dot)/d T
        local_Jac
            .template block<C_size, temperature_size>(C_index,
                                                      temperature_index)
            .noalias() += NpT * ip_cv.dfC_4_MCu_dT * div_u_dot * NT * w;

        GlobalDimMatrixType const advection_C_G = ip.rhoCGR * k_over_mu_G;
        GlobalDimMatrixType const advection_C_L = ip.rhoCLR * k_over_mu_L;
        GlobalDimMatrixType const diffusion_C_G_p =
            -phi_G * ip.rhoGR * D_C_G * ip.dxmWG_dpGR;
        GlobalDimMatrixType const diffusion_C_L_p =
            -phi_L * ip.rhoLR * D_C_L * ip.dxmWL_dpLR;
        GlobalDimMatrixType const diffusion_C_G_T =
            -phi_G * ip.rhoGR * D_C_G * ip.dxmWG_dT;
        GlobalDimMatrixType const diffusion_C_L_T =
            -phi_L * ip.rhoLR * D_C_L * ip.dxmWL_dT;

        GlobalDimMatrixType const advection_C = advection_C_G + advection_C_L;
        GlobalDimMatrixType const diffusion_C_p =
            diffusion_C_G_p + diffusion_C_L_p;
        GlobalDimMatrixType const diffusion_C_T =
            diffusion_C_G_T + diffusion_C_L_T;

        LCpG.noalias() += gradNpT * (advection_C + diffusion_C_p) * gradNp * w;

        // d (fC_4_LCpG * grad p_GR)/d p_GR
        local_Jac.template block<C_size, C_size>(C_index, C_index).noalias() +=
            gradNpT *
            (ip_cv.dadvection_C_dp_GR
             // + ip_cv.ddiffusion_C_p_dp_GR TODO (naumov)
             ) *
            gradpGR * Np * w;

        // d (fC_4_LCpG * grad p_GR)/d p_cap
        local_Jac.template block<C_size, W_size>(C_index, W_index).noalias() +=
            gradNpT *
            (ip_cv.dadvection_C_dp_cap
             // + ip_cv.ddiffusion_C_p_dp_GR TODO (naumov)
             ) *
            gradpGR * Np * w;

        // d (fC_4_LCpG * grad p_GR)/d T
        local_Jac
            .template block<C_size, temperature_size>(C_index,
                                                      temperature_index)
            .noalias() += gradNpT * ip_cv.dfC_4_LCpG_dT * gradpGR * NT * w;

        // d (fC_4_MCpG * p_GR_dot)/d p_GR
        local_Jac.template block<C_size, C_size>(C_index, C_index).noalias() +=
            NpT * ip_cv.dfC_4_MCpG_dp_GR * (pGR - pGR_prev) / dt * Np * w;

        // d (fC_4_MCpG * p_GR_dot)/d T
        local_Jac
            .template block<C_size, temperature_size>(C_index,
                                                      temperature_index)
            .noalias() +=
            NpT * ip_cv.dfC_4_MCpG_dT * (pGR - pGR_prev) / dt * NT * w;

        LCpC.noalias() -=
            gradNpT * (advection_C_L + diffusion_C_L_p) * gradNp * w;

        /* TODO (naumov) This part is not tested by any of the current ctests.
        // d (fC_4_LCpC * grad p_cap)/d p_GR
        local_Jac.template block<C_size, C_size>(C_index, C_index).noalias() +=
            gradNpT *
            (ip_cv.dfC_4_LCpC_a_dp_GR
             // + ip_cv.dfC_4_LCpC_d_dp_GR TODO (naumov)
             ) *
            gradpCap * Np * w;
        // d (fC_4_LCpC * grad p_cap)/d p_cap
        local_Jac.template block<C_size, W_size>(C_index, W_index).noalias() +=
            gradNpT *
            (ip_cv.dfC_4_LCpC_a_dp_cap
             // + ip_cv.dfC_4_LCpC_d_dp_cap TODO (naumov)
             ) *
            gradpCap * Np * w;

        local_Jac
            .template block<C_size, temperature_size>(C_index,
                                                      temperature_index)
            .noalias() += gradNpT *
                          (ip_cv.dfC_4_LCpC_a_dT
                           // + ip_cv.dfC_4_LCpC_d_dT TODO (naumov)
                           ) *
                          gradpCap * Np * w;
        */

        LCT.noalias() += gradNpT * diffusion_C_T * gradNp * w;

        // fC_1
        fC.noalias() += gradNpT *
                        (advection_C_G * ip.rhoGR + advection_C_L * ip.rhoLR) *
                        b * w;

        if (!this->process_data_.apply_mass_lumping)
        {
            // fC_2 = \int a * s_L_dot
            auto const a = phi * (ip.rhoCLR - ip.rhoCGR) -
                           rho_C_FR * pCap * (alpha_B - phi) * beta_p_SR;
            fC.noalias() -= NpT * a * s_L_dot * w;

            local_Jac.template block<C_size, C_size>(C_index, C_index)
                .noalias() +=
                NpT *
                (ip_cv.dfC_2a_dp_GR * s_L_dot /*- a * (ds_L_dp_GR = 0) / dt*/) *
                Np * w;

            local_Jac.template block<C_size, W_size>(C_index, W_index)
                .noalias() +=
                NpT *
                (ip_cv.dfC_2a_dp_cap * s_L_dot + a * ip_cv.dS_L_dp_cap() / dt) *
                Np * w;

            local_Jac
                .template block<C_size, temperature_size>(C_index,
                                                          temperature_index)
                .noalias() += NpT * ip_cv.dfC_2a_dT * s_L_dot * NT * w;
        }
        {
            // fC_3 = \int phi * a
            double const a = s_G * rho_C_GR_dot + s_L * rho_C_LR_dot;
            fC.noalias() -= NpT * phi * a * w;

            local_Jac.template block<C_size, C_size>(C_index, C_index)
                .noalias() += NpT * phi * ip_cv.dfC_3a_dp_GR * Np * w;

            local_Jac.template block<C_size, W_size>(C_index, W_index)
                .noalias() += NpT * phi * ip_cv.dfC_3a_dp_cap * Np * w;

            local_Jac
                .template block<C_size, temperature_size>(C_index,
                                                          temperature_index)
                .noalias() += NpT *
                              (
#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
                                  ip.dphi_dT * a +
#endif  // NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
                                  phi * ip_cv.dfC_3a_dT) *
                              NT * w;
        }
        // ---------------------------------------------------------------------
        // W-component equation
        // ---------------------------------------------------------------------

        MWpG.noalias() += NpT * rho_W_FR * (alpha_B - phi) * beta_p_SR * Np * w;
        MWpC.noalias() -=
            NpT * rho_W_FR * (alpha_B - phi) * beta_p_SR * s_L * Np * w;

        if (this->process_data_.apply_mass_lumping)
        {
            if (pCap - pCap_prev != 0.)  // avoid division by Zero
            {
                MWpC.noalias() +=
                    NpT *
                    (phi * (ip.rhoWLR - ip.rhoWGR) -
                     rho_W_FR * pCap * (alpha_B - phi) * beta_p_SR) *
                    s_L_dot * dt / (pCap - pCap_prev) * Np * w;
            }
        }

        MWT.noalias() -= NpT * rho_W_FR * (alpha_B - phi) * beta_T_SR * Np * w;

        MWu.noalias() += NpT * rho_W_FR * alpha_B * mT * Bu * w;

        GlobalDimMatrixType const advection_W_G = ip.rhoWGR * k_over_mu_G;
        GlobalDimMatrixType const advection_W_L = ip.rhoWLR * k_over_mu_L;
        GlobalDimMatrixType const diffusion_W_G_p =
            phi_G * ip.rhoGR * D_W_G * ip.dxmWG_dpGR;
        GlobalDimMatrixType const diffusion_W_L_p =
            phi_L * ip.rhoLR * D_W_L * ip.dxmWL_dpLR;
        GlobalDimMatrixType const diffusion_W_G_T =
            phi_G * ip.rhoGR * D_W_G * ip.dxmWG_dT;
        GlobalDimMatrixType const diffusion_W_L_T =
            phi_L * ip.rhoLR * D_W_L * ip.dxmWL_dT;

        GlobalDimMatrixType const advection_W = advection_W_G + advection_W_L;
        GlobalDimMatrixType const diffusion_W_p =
            diffusion_W_G_p + diffusion_W_L_p;
        GlobalDimMatrixType const diffusion_W_T =
            diffusion_W_G_T + diffusion_W_L_T;

        LWpG.noalias() += gradNpT * (advection_W + diffusion_W_p) * gradNp * w;

        // fW_4 LWpG' parts; LWpG = \int grad (a + d) grad
        local_Jac.template block<W_size, C_size>(W_index, C_index).noalias() +=
            gradNpT * (ip_cv.dfW_4_LWpG_a_dp_GR + ip_cv.dfW_4_LWpG_d_dp_GR) *
            gradpGR * Np * w;

        local_Jac.template block<W_size, W_size>(W_index, W_index).noalias() +=
            gradNpT * (ip_cv.dfW_4_LWpG_a_dp_cap + ip_cv.dfW_4_LWpG_d_dp_cap) *
            gradpGR * Np * w;

        local_Jac
            .template block<W_size, temperature_size>(W_index,
                                                      temperature_index)
            .noalias() += gradNpT *
                          (ip_cv.dfW_4_LWpG_a_dT + ip_cv.dfW_4_LWpG_d_dT) *
                          gradpGR * NT * w;

        LWpC.noalias() -=
            gradNpT * (advection_W_L + diffusion_W_L_p) * gradNp * w;

        // fW_4 LWp_cap' parts; LWpC = \int grad (a + d) grad
        local_Jac.template block<W_size, C_size>(W_index, C_index).noalias() -=
            gradNpT * (ip_cv.dfW_4_LWpC_a_dp_GR + ip_cv.dfW_4_LWpC_d_dp_GR) *
            gradpCap * Np * w;

        local_Jac.template block<W_size, W_size>(W_index, W_index).noalias() -=
            gradNpT * (ip_cv.dfW_4_LWpC_a_dp_cap + ip_cv.dfW_4_LWpC_d_dp_cap) *
            gradpCap * Np * w;

        local_Jac
            .template block<W_size, temperature_size>(W_index,
                                                      temperature_index)
            .noalias() -= gradNpT *
                          (ip_cv.dfW_4_LWpC_a_dT + ip_cv.dfW_4_LWpC_d_dT) *
                          gradpCap * NT * w;

        LWT.noalias() += gradNpT * (diffusion_W_T)*gradNp * w;

        // fW_1
        fW.noalias() += gradNpT *
                        (advection_W_G * ip.rhoGR + advection_W_L * ip.rhoLR) *
                        b * w;

        // fW_2 = \int (f - g) * s_L_dot
        if (!this->process_data_.apply_mass_lumping)
        {
            double const f = phi * (ip.rhoWLR - ip.rhoWGR);
            double const g = rho_W_FR * pCap * (alpha_B - phi) * beta_p_SR;

            fW.noalias() -= NpT * (f - g) * s_L_dot * w;

            local_Jac.template block<W_size, C_size>(W_index, C_index)
                .noalias() += NpT * (ip_cv.dfW_2a_dp_GR - ip_cv.dfW_2b_dp_GR) *
                              s_L_dot * Np * w;

            // sign negated because of dp_cap = -dp_LR
            // TODO (naumov) Had to change the sign to get equal Jacobian WW
            // blocks in A2 Test. Where is the error?
            local_Jac.template block<W_size, W_size>(W_index, W_index)
                .noalias() +=
                NpT *
                ((ip_cv.dfW_2a_dp_cap - ip_cv.dfW_2b_dp_cap) * s_L_dot +
                 (f - g) * ip_cv.dS_L_dp_cap() / dt) *
                Np * w;

            local_Jac
                .template block<W_size, temperature_size>(W_index,
                                                          temperature_index)
                .noalias() +=
                NpT * (ip_cv.dfW_2a_dT - ip_cv.dfW_2b_dT) * s_L_dot * Np * w;
        }

        // fW_3 = \int phi * a
        fW.noalias() -=
            NpT * phi * (s_G * rho_W_GR_dot + s_L * rho_W_LR_dot) * w;

        local_Jac.template block<W_size, C_size>(W_index, C_index).noalias() +=
            NpT * phi * ip_cv.dfW_3a_dp_GR * Np * w;

        local_Jac.template block<W_size, W_size>(W_index, W_index).noalias() +=
            NpT * phi * ip_cv.dfW_3a_dp_cap * Np * w;

        local_Jac
            .template block<W_size, temperature_size>(W_index,
                                                      temperature_index)
            .noalias() +=
            NpT *
            (
#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
                ip.dphi_dT * (s_G * rho_W_GR_dot + s_L * rho_W_LR_dot) +
#endif  // NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
                phi * ip_cv.dfW_3a_dT) *
            NT * w;

        // ---------------------------------------------------------------------
        //  - temperature equation
        // ---------------------------------------------------------------------

        MTu.noalias() += NTT * rho_h_eff * mT * Bu * w;

        // dfT_4/dp_GR
        // d (MTu * u_dot)/dp_GR
        local_Jac
            .template block<temperature_size, C_size>(temperature_index,
                                                      C_index)
            .noalias() += NTT * ip_cv.drho_h_eff_dp_GR * div_u_dot * NT * w;

        // dfT_4/dp_cap
        // d (MTu * u_dot)/dp_cap
        local_Jac
            .template block<temperature_size, W_size>(temperature_index,
                                                      W_index)
            .noalias() -= NTT * ip_cv.drho_h_eff_dp_cap * div_u_dot * NT * w;

        // dfT_4/dT
        // d (MTu * u_dot)/dT
        local_Jac
            .template block<temperature_size, temperature_size>(
                temperature_index, temperature_index)
            .noalias() += NTT * ip_cv.drho_h_eff_dT * div_u_dot * NT * w;

        KTT.noalias() += gradNTT * ip.lambda * gradNT * w;

        // d KTT/dp_GR * T
        // TODO (naumov) always zero if lambda_xR have no derivatives wrt. p_GR.
        // dlambda_dp_GR =
        //      (dphi_G_dp_GR = 0) * lambdaGR + phi_G * dlambda_GR_dp_GR +
        //      (dphi_L_dp_GR = 0) * lambdaLR + phi_L * dlambda_LR_dp_GR +
        //      (dphi_S_dp_GR = 0) * lambdaSR + phi_S * dlambda_SR_dp_GR +
        //      = 0
        //
        // Since dlambda/dp_GR is 0 the derivative is omitted:
        // local_Jac
        //    .template block<temperature_size, C_size>(temperature_index,
        //                                              C_index)
        //    .noalias() += gradNTT * dlambda_dp_GR * gradT * Np * w;

        // d KTT/dp_cap * T
        local_Jac
            .template block<temperature_size, W_size>(temperature_index,
                                                      W_index)
            .noalias() += gradNTT * ip_cv.dlambda_dp_cap * gradT * Np * w;

        // d KTT/dT * T
        local_Jac
            .template block<temperature_size, temperature_size>(
                temperature_index, temperature_index)
            .noalias() += gradNTT * ip_cv.dlambda_dT * gradT * NT * w;

        // fT_1
        fT.noalias() -= NTT * rho_u_eff_dot * w;

        // dfT_1/dp_GR
        local_Jac
            .template block<temperature_size, C_size>(temperature_index,
                                                      C_index)
            .noalias() += NpT / dt * ip_cv.drho_u_eff_dp_GR * Np * w;

        // dfT_1/dp_cap
        local_Jac
            .template block<temperature_size, W_size>(temperature_index,
                                                      W_index)
            .noalias() += NpT / dt * ip_cv.drho_u_eff_dp_cap * Np * w;

        // dfT_1/dT
        // MTT
        local_Jac
            .template block<temperature_size, temperature_size>(
                temperature_index, temperature_index)
            .noalias() += NTT * ip_cv.drho_u_eff_dT / dt * NT * w;

        // fT_2
        fT.noalias() +=
            gradNTT * (ip.rhoGR * h_G * ip.w_GS + ip.rhoLR * h_L * ip.w_LS) * w;

        // dfT_2/dp_GR
        local_Jac
            .template block<temperature_size, C_size>(temperature_index,
                                                      C_index)
            .noalias() -=
            // dfT_2/dp_GR first part
            gradNTT * ip_cv.drho_GR_h_w_eff_dp_GR_Npart * Np * w +
            // dfT_2/dp_GR second part
            gradNTT * ip_cv.drho_GR_h_w_eff_dp_GR_gradNpart * gradNp * w;

        // dfT_2/dp_cap
        local_Jac
            .template block<temperature_size, W_size>(temperature_index,
                                                      W_index)
            .noalias() -=
            // first part of dfT_2/dp_cap
            gradNTT * (-ip_cv.drho_LR_h_w_eff_dp_cap_Npart) * Np * w +
            // second part of dfT_2/dp_cap
            gradNTT * (-ip_cv.drho_LR_h_w_eff_dp_cap_gradNpart) * gradNp * w;

        // dfT_2/dT
        local_Jac
            .template block<temperature_size, temperature_size>(
                temperature_index, temperature_index)
            .noalias() -= gradNTT * ip_cv.drho_GR_h_w_eff_dT * NT * w;

        // fT_3
        fT.noalias() +=
            NTT *
            (ip.rhoGR * ip.w_GS.transpose() + ip.rhoLR * ip.w_LS.transpose()) *
            b * w;

        fT.noalias() +=
            gradNTT *
            (ip.rhoCGR * ip.h_CG * ip.d_CG + ip.rhoWGR * ip.h_WG * ip.d_WG) * w;

        // ---------------------------------------------------------------------
        //  - displacement equation
        // ---------------------------------------------------------------------

        KUpG.noalias() -= (BuT * alpha_B * m * Np) * w;

        // dfU_2/dp_GR = dKUpG/dp_GR * p_GR + KUpG. The former is zero, the
        // latter is handled below.

        KUpC.noalias() += (BuT * alpha_B * ip_cv.chi_S_L.chi_S_L * m * Np) * w;

        // dfU_2/dp_cap = dKUpC/dp_cap * p_cap + KUpC. The former is handled
        // here, the latter below.
        local_Jac
            .template block<displacement_size, W_size>(displacement_index,
                                                       W_index)
            .noalias() += BuT * alpha_B * ip_cv.chi_S_L.dchi_dS_L *
                          ip_cv.dS_L_dp_cap() * pCap * m * Np * w;

        local_Jac
            .template block<displacement_size, displacement_size>(
                displacement_index, displacement_index)
            .noalias() += BuT * ip_cd.s_mech_data.stiffness_tensor * Bu * w;

        // fU_1
        fU.noalias() -= (BuT * current_state.eff_stress_data.sigma -
                         N_u_op(Nu).transpose() * rho * b) *
                        w;

        // KuT
        local_Jac
            .template block<displacement_size, temperature_size>(
                displacement_index, temperature_index)
            .noalias() -=
            BuT *
            (ip_cd.s_mech_data.stiffness_tensor *
             ip_cv.s_therm_exp_data.solid_linear_thermal_expansivity_vector) *
            NT * w;

        /* TODO (naumov) Test with gravity needed to check this Jacobian part.
        local_Jac
            .template block<displacement_size, temperature_size>(
                displacement_index, temperature_index)
            .noalias() += N_u_op(Nu).transpose() * ip_cv.drho_dT * b *
                          N_u_op(Nu).transpose() * w;
         */

        if (this->process_data_.apply_mass_lumping)
        {
            MCpG = MCpG.colwise().sum().eval().asDiagonal();
            MCpC = MCpC.colwise().sum().eval().asDiagonal();
            MWpG = MWpG.colwise().sum().eval().asDiagonal();
            MWpC = MWpC.colwise().sum().eval().asDiagonal();
        }
    }  // int_point-loop

    // --- Gas ---
    // fC_4
    fC.noalias() -= LCpG * gas_pressure + LCpC * capillary_pressure +
                    LCT * temperature +
                    MCpG * (gas_pressure - gas_pressure_prev) / dt +
                    MCpC * (capillary_pressure - capillary_pressure_prev) / dt +
                    MCT * (temperature - temperature_prev) / dt +
                    MCu * (displacement - displacement_prev) / dt;

    local_Jac.template block<C_size, C_size>(C_index, C_index).noalias() +=
        LCpG + MCpG / dt;
    local_Jac.template block<C_size, W_size>(C_index, W_index).noalias() +=
        LCpC + MCpC / dt;
    local_Jac
        .template block<C_size, temperature_size>(C_index, temperature_index)
        .noalias() += LCT + MCT / dt;
    local_Jac
        .template block<C_size, displacement_size>(C_index, displacement_index)
        .noalias() += MCu / dt;

    // --- Capillary pressure ---
    // fW_4
    fW.noalias() -= LWpG * gas_pressure + LWpC * capillary_pressure +
                    LWT * temperature +
                    MWpG * (gas_pressure - gas_pressure_prev) / dt +
                    MWpC * (capillary_pressure - capillary_pressure_prev) / dt +
                    MWT * (temperature - temperature_prev) / dt +
                    MWu * (displacement - displacement_prev) / dt;

    local_Jac.template block<W_size, W_size>(W_index, W_index).noalias() +=
        LWpC + MWpC / dt;
    local_Jac.template block<W_size, C_size>(W_index, C_index).noalias() +=
        LWpG + MWpG / dt;
    local_Jac
        .template block<W_size, temperature_size>(W_index, temperature_index)
        .noalias() += LWT + MWT / dt;
    local_Jac
        .template block<W_size, displacement_size>(W_index, displacement_index)
        .noalias() += MWu / dt;

    // --- Temperature ---
    // fT_4
    fT.noalias() -=
        KTT * temperature + MTu * (displacement - displacement_prev) / dt;

    local_Jac
        .template block<temperature_size, temperature_size>(temperature_index,
                                                            temperature_index)
        .noalias() += KTT;
    local_Jac
        .template block<temperature_size, displacement_size>(temperature_index,
                                                             displacement_index)
        .noalias() += MTu / dt;

    // --- Displacement ---
    // fU_2
    fU.noalias() -= KUpG * gas_pressure + KUpC * capillary_pressure;

    local_Jac
        .template block<displacement_size, C_size>(displacement_index, C_index)
        .noalias() += KUpG;
    local_Jac
        .template block<displacement_size, W_size>(displacement_index, W_index)
        .noalias() += KUpC;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> const& TH2MLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    getIntPtDarcyVelocityGas(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();

    cache.clear();
    auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, DisplacementDim, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        cache_matrix.col(ip).noalias() = _ip_data[ip].w_GS;
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> const& TH2MLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    getIntPtDarcyVelocityLiquid(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();

    cache.clear();
    auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, DisplacementDim, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        cache_matrix.col(ip).noalias() = _ip_data[ip].w_LS;
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> const& TH2MLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    getIntPtDiffusionVelocityVapourGas(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();

    cache.clear();
    auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, DisplacementDim, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        cache_matrix.col(ip).noalias() = _ip_data[ip].d_WG;
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> const& TH2MLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    getIntPtDiffusionVelocityGasGas(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();

    cache.clear();
    auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, DisplacementDim, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        cache_matrix.col(ip).noalias() = _ip_data[ip].d_CG;
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> const& TH2MLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    getIntPtDiffusionVelocitySoluteLiquid(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();

    cache.clear();
    auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, DisplacementDim, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        cache_matrix.col(ip).noalias() = _ip_data[ip].d_CL;
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> const& TH2MLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    getIntPtDiffusionVelocityLiquidLiquid(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();

    cache.clear();
    auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, DisplacementDim, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        cache_matrix.col(ip).noalias() = _ip_data[ip].d_WL;
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void TH2MLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                        DisplacementDim>::
    computeSecondaryVariableConcrete(double const t, double const dt,
                                     Eigen::VectorXd const& local_x,
                                     Eigen::VectorXd const& local_x_prev)
{
    auto const gas_pressure =
        local_x.template segment<gas_pressure_size>(gas_pressure_index);
    auto const capillary_pressure =
        local_x.template segment<capillary_pressure_size>(
            capillary_pressure_index);
    auto const liquid_pressure = (gas_pressure - capillary_pressure).eval();

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(this->element_, this->is_axially_symmetric_,
                         gas_pressure,
                         *this->process_data_.gas_pressure_interpolated);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(this->element_, this->is_axially_symmetric_,
                         capillary_pressure,
                         *this->process_data_.capillary_pressure_interpolated);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(this->element_, this->is_axially_symmetric_,
                         liquid_pressure,
                         *this->process_data_.liquid_pressure_interpolated);

    auto const temperature =
        local_x.template segment<temperature_size>(temperature_index);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(this->element_, this->is_axially_symmetric_,
                         temperature,
                         *this->process_data_.temperature_interpolated);

    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();

    double saturation_avg = 0;

    ConstitutiveRelations::ConstitutiveModels<DisplacementDim> const models{
        this->solid_material_};

    updateConstitutiveVariables(local_x, local_x_prev, t, dt, models);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        saturation_avg += this->current_states_[ip].S_L_data.S_L;
    }
    saturation_avg /= n_integration_points;
    (*this->process_data_.element_saturation)[this->element_.getID()] =
        saturation_avg;
}

}  // namespace TH2M
}  // namespace ProcessLib
