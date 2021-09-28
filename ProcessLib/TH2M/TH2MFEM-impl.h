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

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"

namespace ProcessLib
{
namespace TH2M
{
namespace MPL = MaterialPropertyLib;

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
TH2MLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                   IntegrationMethod, DisplacementDim>::
    TH2MLocalAssembler(MeshLib::Element const& e,
                       std::size_t const /*local_matrix_size*/,
                       bool const is_axially_symmetric,
                       unsigned const integration_order,
                       TH2MProcessData<DisplacementDim>& process_data)
    : _process_data(process_data),
      _integration_method(integration_order),
      _element(e),
      _is_axially_symmetric(is_axially_symmetric)
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);
    _secondary_data.N_u.resize(n_integration_points);

    auto const shape_matrices_u =
        NumLib::initShapeMatrices<ShapeFunctionDisplacement,
                                  ShapeMatricesTypeDisplacement,
                                  DisplacementDim>(e, is_axially_symmetric,
                                                   _integration_method);

    auto const shape_matrices_p =
        NumLib::initShapeMatrices<ShapeFunctionPressure,
                                  ShapeMatricesTypePressure, DisplacementDim>(
            e, is_axially_symmetric, _integration_method);

    auto const& solid_material =
        MaterialLib::Solids::selectSolidConstitutiveRelation(
            _process_data.solid_materials, _process_data.material_ids,
            e.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        _ip_data.emplace_back(solid_material);
        auto& ip_data = _ip_data[ip];
        auto const& sm_u = shape_matrices_u[ip];
        ip_data.integration_weight =
            _integration_method.getWeightedPoint(ip).getWeight() *
            sm_u.integralMeasure * sm_u.detJ;

        ip_data.N_u_op = ShapeMatricesTypeDisplacement::template MatrixType<
            DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                      displacement_size);
        for (int i = 0; i < DisplacementDim; ++i)
        {
            ip_data.N_u_op
                .template block<1, displacement_size / DisplacementDim>(
                    i, i * displacement_size / DisplacementDim)
                .noalias() = sm_u.N;
        }

        ip_data.N_u = sm_u.N;
        ip_data.dNdx_u = sm_u.dNdx;

        ip_data.N_p = shape_matrices_p[ip].N;
        ip_data.dNdx_p = shape_matrices_p[ip].dNdx;

        _secondary_data.N_u[ip] = shape_matrices_u[ip].N;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<ConstitutiveVariables<DisplacementDim>>
TH2MLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                   IntegrationMethod, DisplacementDim>::
    updateConstitutiveVariables(Eigen::VectorXd const& local_x,
                                Eigen::VectorXd const& local_x_dot,
                                double const t, double const dt)
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
    auto const temperature_dot =
        local_x_dot.template segment<temperature_size>(temperature_index);

    auto const displacement =
        local_x.template segment<displacement_size>(displacement_index);

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());

    auto const& medium = *_process_data.media_map->getMedium(_element.getID());
    auto const& gas_phase = medium.phase("Gas");
    auto const& liquid_phase = medium.phase("AqueousLiquid");
    auto const& solid_phase = medium.phase("Solid");

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    std::vector<ConstitutiveVariables<DisplacementDim>>
        ip_constitutive_variables(n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);
        auto& ip_data = _ip_data[ip];
        auto& ip_cv = ip_constitutive_variables[ip];

        auto const& Np = ip_data.N_p;
        auto const& NT = Np;
        auto const& Nu = ip_data.N_u;
        auto const& gradNu = ip_data.dNdx_u;
        auto const& gradNp = ip_data.dNdx_p;
        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                _element, Nu);

        double const T = NT.dot(temperature);
        double const pGR = Np.dot(gas_pressure);
        double const pCap = Np.dot(capillary_pressure);
        double const pLR = pGR - pCap;
        GlobalDimVectorType const gradpGR = gradNp * gas_pressure;
        GlobalDimVectorType const gradpCap = gradNp * capillary_pressure;

        MPL::VariableArray vars;
        vars[static_cast<int>(MPL::Variable::temperature)] = T;
        vars[static_cast<int>(MPL::Variable::phase_pressure)] = pGR;
        vars[static_cast<int>(MPL::Variable::capillary_pressure)] = pCap;
        vars[static_cast<int>(MPL::Variable::liquid_phase_pressure)] = pLR;

        // medium properties
        auto const K_S = ip_data.solid_material.getBulkModulus(t, pos);

        ip_data.alpha_B = medium.property(MPL::PropertyType::biot_coefficient)
                              .template value<double>(vars, pos, t, dt);

        ip_data.s_L =
            medium.property(MPL::PropertyType::saturation)
                .template value<double>(
                    vars, pos, t, std::numeric_limits<double>::quiet_NaN());

        vars[static_cast<int>(MPL::Variable::liquid_saturation)] = ip_data.s_L;

        // intrinsic permeability
        ip_data.k_S = MPL::formEigenTensor<DisplacementDim>(
            medium.property(MPL::PropertyType::permeability)
                .value(vars, pos, t, dt));

        auto const Bu =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                gradNu, Nu, x_coord, _is_axially_symmetric);

        auto& eps = ip_data.eps;
        eps.noalias() = Bu * displacement;

        // relative permeability
        // Set mechanical variables for the intrinsic permeability model
        // For stress dependent permeability.
        {
            // Note: if Bishop model is available, ip_data.s_L in the following
            // computation should be replaced with the Bishop value.
            auto const sigma_total =
                (_ip_data[ip].sigma_eff - ip_data.alpha_B *
                                              (pGR - ip_data.s_L * pCap) *
                                              Invariants::identity2)
                    .eval();

            vars[static_cast<int>(MPL::Variable::total_stress)]
                .emplace<SymmetricTensor>(
                    MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                        sigma_total));
        }
        // For strain dependent permeability
        vars[static_cast<int>(MPL::Variable::volumetric_strain)] =
            Invariants::trace(eps);
        vars[static_cast<int>(MPL::Variable::equivalent_plastic_strain)] =
            _ip_data[ip].material_state_variables->getEquivalentPlasticStrain();
        vars[static_cast<int>(MPL::Variable::mechanical_strain)]
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps);

        ip_data.k_rel_G =
            medium
                .property(
                    MPL::PropertyType::relative_permeability_nonwetting_phase)
                .template value<double>(vars, pos, t, dt);

        ip_data.k_rel_L =
            medium.property(MPL::PropertyType::relative_permeability)
                .template value<double>(vars, pos, t, dt);

        // solid phase compressibility
        ip_data.beta_p_SR = (1. - ip_data.alpha_B) / K_S;

        // solid phase linear thermal expansion coefficient
        ip_data.alpha_T_SR = MathLib::KelvinVector::tensorToKelvin<
            DisplacementDim>(MaterialPropertyLib::formEigenTensor<3>(
            solid_phase
                .property(
                    MaterialPropertyLib::PropertyType::thermal_expansivity)
                .value(vars, pos, t, dt)));

        // isotropic solid phase volumetric thermal expansion coefficient
        ip_data.beta_T_SR = Invariants::trace(ip_data.alpha_T_SR);

        double const T_dot = NT.dot(temperature_dot);
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const
            dthermal_strain = ip_data.alpha_T_SR * T_dot * dt;

        auto& eps_prev = ip_data.eps_prev;
        auto& eps_m = ip_data.eps_m;
        auto& eps_m_prev = ip_data.eps_m_prev;

        eps_m.noalias() = eps_m_prev + eps - eps_prev - dthermal_strain;
        vars[static_cast<int>(MaterialPropertyLib::Variable::mechanical_strain)]
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m);

        auto const rho_ref_SR =
            solid_phase.property(MPL::PropertyType::density)
                .template value<double>(
                    vars, pos, t, std::numeric_limits<double>::quiet_NaN());

        auto const lambdaSR = MPL::formEigenTensor<DisplacementDim>(
            solid_phase.property(MPL::PropertyType::thermal_conductivity)
                .value(vars, pos, t, dt));

        double const T0 = _process_data.reference_temperature(t, pos)[0];
        double const delta_T(T - T0);
        ip_data.thermal_volume_strain = ip_data.beta_T_SR * delta_T;

        // initial porosity
        auto const phi_0 = medium.property(MPL::PropertyType::porosity)
                               .template value<double>(vars, pos, t, dt);

        auto const phi_S_0 = 1. - phi_0;

#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
        auto const& m = Invariants::identity2;
        double const div_u = m.transpose() * eps;

        const double phi_S = phi_S_0 * (1. + ip_data.thermal_volume_strain -
                                        ip_data.alpha_B * div_u);
#else   // NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
        const double phi_S = phi_S_0;
#endif  // NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION

        // porosity
        ip_data.phi = 1. - phi_S;

        // solid phase density
#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
        auto const rhoSR = rho_ref_SR * (1. - ip_data.thermal_volume_strain +
                                         (ip_data.alpha_B - 1.) * div_u);
#else   // NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
        auto const rhoSR = rho_ref_SR;
#endif  // NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION

        auto const T_prev = T - T_dot * dt;
        ip_cv.C = ip_data.updateConstitutiveRelation(vars, t, pos, dt, T_prev);

        // constitutive model object as specified in process creation
        auto& ptm = *_process_data.phase_transition_model_;
        ptm.computeConstitutiveVariables(&medium, vars, pos, t, dt);
        auto& c = ptm.cv;
        auto const phi_L = ip_data.s_L * ip_data.phi;
        auto const phi_G = (1. - ip_data.s_L) * ip_data.phi;

        // TODO (Grunwald): individual volume fractions can be stored in a
        // container of type MPL::Composition (a.k.a. std::vector<double>) which
        // can be stored in the variable array for access in MPL properties
        // ---
        // MaterialPropertyLib::Composition volume_fraction{phi_G, phi_L,
        // phi_S};
        // vars[static_cast<int>(MPL::Variable::volume_fraction)] =
        //     volume_fraction;
        // ---
        // TODO (Grunwald) replace effective thermal conductivity by a more
        // sophisticated law by allowing the law to be chosen in the project
        // file as medium property, e.g.
        // lambda = medium.property(MPL::PropertyType::thermal_conductivity)..
        // where volume fraction is stored in the variable array

        auto const lambdaGR = MPL::formEigenTensor<DisplacementDim>(c.lambdaGR);
        auto const lambdaLR = MPL::formEigenTensor<DisplacementDim>(c.lambdaLR);

        ip_data.lambda = phi_G * lambdaGR + phi_L * lambdaLR + phi_S * lambdaSR;

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

        ip_data.dxmCG_dpGR = c.dxmCG_dpGR;
        ip_data.dxmCG_dT = c.dxmCG_dT;
        ip_data.dxmCL_dpLR = c.dxmCL_dpLR;
        ip_data.dxmCL_dT = c.dxmCL_dT;
        ip_data.dxmWG_dpGR = c.dxmWG_dpGR;
        ip_data.dxmWG_dT = c.dxmWG_dT;
        ip_data.dxmWL_dpLR = c.dxmWL_dpLR;
        ip_data.dxmWL_dT = c.dxmWL_dT;

        // for variable output
        ip_data.xnCG = c.xnCG;
        ip_data.xmCG = c.xmCG;
        ip_data.xmWL = c.xmWL;

        ip_data.diffusion_coefficient_vapour = c.diffusion_coefficient_vapour;
        ip_data.diffusion_coefficient_solvate = c.diffusion_coefficient_solvate;

        ip_data.h_G = c.hG;
        ip_data.h_L = c.hL;
        ip_data.pWGR = c.pWGR;

        // ---------------------------------------------------------------------
        // Derivatives for Jacobian
        // ---------------------------------------------------------------------
        ip_cv.drho_LR_dT =
            liquid_phase.property(MPL::PropertyType::density)
                .template dValue<double>(vars, MPL::Variable::temperature, pos,
                                         t, dt);
        ip_cv.drho_SR_dT =
            solid_phase.property(MPL::PropertyType::density)
                    .template dValue<double>(vars, MPL::Variable::temperature,
                                             pos, t, dt)
#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
                * (1. - ip_data.thermal_volume_strain +
                   (ip_data.alpha_B - 1.) * div_u) -
            rho_ref_SR * ip_data.beta_T_SR
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
                                        ip_data.alpha_B * div_u) +
                                 phi_S_0 * ip_data.beta_T_SR
#endif
            ;

        ip_cv.drho_u_eff_dT =
            phi_G * c.drho_GR_dT * c.uG + phi_G * c.rhoGR * c.du_G_dT +
            phi_L * ip_cv.drho_LR_dT * c.uL + phi_L * c.rhoLR * c.du_L_dT +
            phi_S * ip_cv.drho_SR_dT * u_S + phi_S * rhoSR * cpS +
            dphi_S_dT * rhoSR * u_S;

        ip_cv.ds_L_dp_cap =
            medium[MPL::PropertyType::saturation].template dValue<double>(
                vars, MPL::Variable::capillary_pressure, pos, t, dt);
        // TODO (naumov) Extend for partially saturated media.
        constexpr double ds_L_dp_GR = 0;
        double const ds_G_dp_GR = -ds_L_dp_GR;
        double const ds_G_dp_cap = -ip_cv.ds_L_dp_cap;

        double const dphi_G_dp_GR = -ds_L_dp_GR * ip_data.phi;
        double const dphi_G_dp_cap = -ip_cv.ds_L_dp_cap * ip_data.phi;
        double const dphi_L_dp_GR = ds_L_dp_GR * ip_data.phi;
        double const dphi_L_dp_cap = ip_cv.ds_L_dp_cap * ip_data.phi;

        auto const dlambda_GR_dT = MPL::formEigenTensor<DisplacementDim>(
            gas_phase[MPL::PropertyType::thermal_conductivity].dValue(
                vars, MPL::Variable::temperature, pos, t, dt));
        auto const dlambda_LR_dT = MPL::formEigenTensor<DisplacementDim>(
            liquid_phase[MPL::PropertyType::thermal_conductivity].dValue(
                vars, MPL::Variable::temperature, pos, t, dt));
        auto const dlambda_SR_dT = MPL::formEigenTensor<DisplacementDim>(
            solid_phase[MPL::PropertyType::thermal_conductivity].dValue(
                vars, MPL::Variable::temperature, pos, t, dt));

        ip_cv.dlambda_dp_GR = dphi_G_dp_GR * lambdaGR + dphi_L_dp_GR * lambdaLR;

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
        // And also for drho_GR/dp_LR = drho_GR/dp_GR * dp_GR/dp_LR:
        double const drho_GR_dp_cap = 0;

        ip_cv.drho_h_eff_dp_GR =
            dphi_G_dp_GR * c.rhoGR * c.hG + phi_G * c.drho_GR_dp_GR * c.hG +
            dphi_L_dp_GR * c.rhoLR * c.hL + phi_L * drho_LR_dp_GR * c.hL;
        ip_cv.drho_h_eff_dp_cap =
            dphi_G_dp_cap * c.rhoGR * c.hG + phi_G * drho_GR_dp_cap * c.hG +
            dphi_L_dp_cap * c.rhoLR * c.hL + phi_L * drho_LR_dp_cap * c.hL;

        // TODO (naumov) Extend for temperature dependent porosities.
        constexpr double dphi_G_dT = 0;
        constexpr double dphi_L_dT = 0;
        ip_cv.drho_h_eff_dT =
            dphi_G_dT * c.rhoGR * c.hG + phi_G * c.drho_GR_dT * c.hG +
            phi_G * c.rhoGR * c.dh_G_dT + dphi_L_dT * c.rhoLR * c.hL +
            phi_L * ip_cv.drho_LR_dT * c.hL + phi_L * c.rhoLR * c.dh_L_dT +
            dphi_S_dT * rhoSR * ip_data.h_S +
            phi_S * ip_cv.drho_SR_dT * ip_data.h_S + phi_S * rhoSR * cpS;

        ip_cv.drho_u_eff_dp_GR =
            dphi_G_dp_GR * c.rhoGR * c.uG + phi_G * c.drho_GR_dp_GR * c.uG +
            phi_G * c.rhoGR * c.du_G_dp_GR + dphi_L_dp_GR * c.rhoLR * c.uL +
            phi_L * drho_LR_dp_GR * c.uL + phi_L * c.rhoLR * c.du_L_dp_GR;

        ip_cv.drho_u_eff_dp_cap =
            dphi_G_dp_cap * c.rhoGR * c.uG + phi_G * drho_GR_dp_cap * c.uG +
            dphi_L_dp_cap * c.rhoLR * c.uL + phi_L * drho_LR_dp_cap * c.uL +
            phi_L * c.rhoLR * c.du_L_dp_cap;

        auto const& b = _process_data.specific_body_force;
        auto const k_over_mu_G =
            (ip_data.k_S * ip_data.k_rel_G / ip_data.muGR).eval();
        auto const k_over_mu_L =
            (ip_data.k_S * ip_data.k_rel_L / ip_data.muLR).eval();
        GlobalDimVectorType const w_GS =
            k_over_mu_G * c.rhoGR * b - k_over_mu_G * gradpGR;
        GlobalDimVectorType const w_LS = k_over_mu_L * gradpCap +
                                         k_over_mu_L * c.rhoLR * b -
                                         k_over_mu_L * gradpGR;

        ip_cv.drho_GR_h_w_eff_dp_GR_Npart =
            c.drho_GR_dp_GR * c.hG * w_GS +
            c.rhoGR * c.hG * k_over_mu_G * c.drho_GR_dp_GR * b;
        ip_cv.drho_GR_h_w_eff_dp_GR_gradNpart =
            -c.rhoGR * c.hG * k_over_mu_G - c.rhoLR * c.hL * k_over_mu_L;

        ip_cv.drho_LR_h_w_eff_dp_cap_Npart =
            -drho_LR_dp_cap * c.hL * w_LS -
            c.rhoLR * c.hL * k_over_mu_L * drho_LR_dp_cap * b;
        ip_cv.drho_LR_h_w_eff_dp_cap_gradNpart =
            // TODO (naumov) why the minus sign??????
            -c.rhoLR * c.hL * k_over_mu_L;

        ip_cv.drho_GR_h_w_eff_dT =
            c.drho_GR_dT * c.hG * w_GS + c.rhoGR * c.dh_G_dT * w_GS +
            ip_cv.drho_LR_dT * c.hL * w_LS + c.rhoLR * c.dh_L_dT * w_LS;
        // TODO (naumov) + k_over_mu_G * drho_GR_dT * b + k_over_mu_L *
        // drho_LR_dT * b

        // Derivatives of s_G * rho_C_GR_dot + s_L * rho_C_LR_dot abbreviated
        // here with S_rho_C_eff.
        double const s_L = ip_data.s_L;
        double const s_G = 1. - ip_data.s_L;
        double const rho_C_GR_dot = (ip_data.rhoCGR - ip_data.rhoCGR_prev) / dt;
        double const rho_C_LR_dot = (ip_data.rhoCLR - ip_data.rhoCLR_prev) / dt;
        double const rho_C_FR = s_G * ip_data.rhoCGR + s_L * ip_data.rhoCLR;
        double const rho_W_FR = s_G * ip_data.rhoWGR + s_L * ip_data.rhoWLR;
        // TODO (naumov) Extend for partially saturated media.
        constexpr double drho_C_GR_dp_cap = 0;
        ip_cv.dfC_3a_dp_GR =
            ds_G_dp_GR * rho_C_GR_dot + s_G * c.drho_C_GR_dp_GR / dt +
            ds_L_dp_GR * rho_C_LR_dot + s_L * c.drho_C_LR_dp_GR / dt;
        ip_cv.dfC_3a_dp_cap =
            ds_G_dp_cap * rho_C_GR_dot + s_G * drho_C_GR_dp_cap / dt +
            ip_cv.ds_L_dp_cap * rho_C_LR_dot - s_L * c.drho_C_LR_dp_LR / dt;
        ip_cv.dfC_3a_dT = s_G * c.drho_C_GR_dT / dt + s_L * c.drho_C_LR_dT / dt;

        double const drho_C_FR_dp_GR =
            ds_G_dp_GR * ip_data.rhoCGR + s_G * c.drho_C_GR_dp_GR +
            ds_L_dp_GR * ip_data.rhoCLR + s_L * c.drho_C_LR_dp_GR;
        ip_cv.dfC_4_MCpG_dp_GR = drho_C_FR_dp_GR *
                                 (ip_data.alpha_B - ip_data.phi) *
                                 ip_data.beta_p_SR;

        double const drho_C_FR_dT = s_G * c.drho_C_GR_dT + s_L * c.drho_C_LR_dT;
        ip_cv.dfC_4_MCpG_dT =
            drho_C_FR_dT * (ip_data.alpha_B - ip_data.phi) * ip_data.beta_p_SR
#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
            - rho_C_FR * ip_data.dphi_dT * ip_data.beta_p_SR
#endif
            ;

        ip_cv.dfC_4_MCT_dT =
            drho_C_FR_dT * (ip_data.alpha_B - ip_data.phi) * ip_data.beta_T_SR
#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
            + rho_C_FR * (ip_data.alpha_B - ip_data.dphi_dT) * ip_data.beta_T_SR
#endif
            ;

        ip_cv.dfC_4_MCu_dT = drho_C_FR_dT * ip_data.alpha_B;

        ip_cv.dfC_2a_dp_GR = -ip_data.phi * c.drho_C_GR_dp_GR -
                             drho_C_FR_dp_GR * pCap *
                                 (ip_data.alpha_B - ip_data.phi) *
                                 ip_data.beta_p_SR;

        double const drho_C_FR_dp_cap =
            ds_G_dp_cap * ip_data.rhoCGR + s_G * drho_C_GR_dp_cap +
            ip_cv.ds_L_dp_cap * ip_data.rhoCLR - s_L * c.drho_C_LR_dp_LR;

        ip_cv.dfC_2a_dp_cap =
            ip_data.phi * (-c.drho_C_LR_dp_LR - drho_C_GR_dp_cap) -
            drho_C_FR_dp_cap * pCap * (ip_data.alpha_B - ip_data.phi) *
                ip_data.beta_p_SR +
            rho_C_FR * (ip_data.alpha_B - ip_data.phi) * ip_data.beta_p_SR;

        ip_cv.dfC_2a_dT =
#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
            ip_data.dphi_dT * (ip_data.rhoCLR - ip_data.rhoCGR) +
#endif
            ip_data.phi * (c.drho_C_LR_dT - c.drho_C_GR_dT) -
            drho_C_FR_dT * pCap * (ip_data.alpha_B - ip_data.phi) *
                ip_data.beta_p_SR
#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
            + rho_C_FR * pCap * ip_data.dphi_dT * ip_data.beta_p_SR
#endif
            ;

        ip_cv.dadvection_C_dp_GR =
            c.drho_C_GR_dp_GR * k_over_mu_G + c.drho_C_LR_dp_GR * k_over_mu_L;

        ip_cv.dfC_4_LCpG_dT =
            c.drho_C_GR_dT * k_over_mu_G + c.drho_C_LR_dT * k_over_mu_L
            // + ip_cv.ddiffusion_C_p_dT TODO (naumov)
            ;

        double const drho_W_FR_dp_GR =
            ds_G_dp_GR * ip_data.rhoWGR + s_G * c.drho_W_GR_dp_GR +
            ds_L_dp_GR * ip_data.rhoWLR + s_L * c.drho_W_LR_dp_GR;
        double const drho_W_FR_dp_cap =
            ds_G_dp_cap * ip_data.rhoWGR + s_G * c.drho_W_GR_dp_cap +
            ip_cv.ds_L_dp_cap * ip_data.rhoWLR - s_L * c.drho_W_LR_dp_LR;
        double const drho_W_FR_dT = s_G * c.drho_W_GR_dT + s_L * c.drho_W_LR_dT;

        ip_cv.dfW_2a_dp_GR =
            ip_data.phi * (c.drho_W_LR_dp_GR - c.drho_W_GR_dp_GR);
        ip_cv.dfW_2b_dp_GR = drho_W_FR_dp_GR * pCap *
                             (ip_data.alpha_B - ip_data.phi) *
                             ip_data.beta_p_SR;
        ip_cv.dfW_2a_dp_cap =
            ip_data.phi * (-c.drho_W_LR_dp_LR - c.drho_W_GR_dp_cap);
        ip_cv.dfW_2b_dp_cap =
            drho_W_FR_dp_cap * pCap * (ip_data.alpha_B - ip_data.phi) *
                ip_data.beta_p_SR +
            rho_W_FR * (ip_data.alpha_B - ip_data.phi) * ip_data.beta_p_SR;

        ip_cv.dfW_2a_dT =
#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
            ip_data.dphi_dT * (ip_data.rhoWLR - ip_data.rhoWGR) +
#endif
            ip_data.phi * (c.drho_W_LR_dT - c.drho_W_GR_dT);
        ip_cv.dfW_2b_dT =
            drho_W_FR_dT * pCap * (ip_data.alpha_B - ip_data.phi) *
                ip_data.beta_p_SR
#ifdef NON_CONSTANT_SOLID_PHASE_VOLUME_FRACTION
            - rho_W_FR * pCap * ip_data.dphi_dT * ip_data.beta_p_SR
#endif
            ;

        double const rho_W_GR_dot = (ip_data.rhoWGR - ip_data.rhoWGR_prev) / dt;
        double const rho_W_LR_dot = (ip_data.rhoWLR - ip_data.rhoWLR_prev) / dt;

        ip_cv.dfW_3a_dp_GR =
            ds_G_dp_GR * rho_W_GR_dot + s_G * c.drho_W_GR_dp_GR / dt +
            ds_L_dp_GR * rho_W_LR_dot + s_L * c.drho_W_LR_dp_GR / dt;
        ip_cv.dfW_3a_dp_cap =
            ds_G_dp_cap * rho_W_GR_dot + s_G * c.drho_W_GR_dp_cap / dt +
            ip_cv.ds_L_dp_cap * rho_W_LR_dot - s_L * c.drho_W_LR_dp_LR / dt;
        ip_cv.dfW_3a_dT = s_G * c.drho_W_GR_dT / dt + s_L * c.drho_W_LR_dT / dt;

        ip_cv.dfW_4a_dp_GR = c.drho_W_GR_dp_GR * k_over_mu_G;
        ip_cv.dfW_4a_dp_cap = -c.drho_W_LR_dp_LR * k_over_mu_L;
        ip_cv.dfW_4a_dT =
            c.drho_W_GR_dT * k_over_mu_G + c.drho_W_LR_dT * k_over_mu_L;

        // TODO (naumov) for dxmW*/d* != 0
        ip_cv.dfW_4d_dp_GR =
            Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Zero();
        ip_cv.dfW_4d_dp_cap =
            Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Zero();
        ip_cv.dfW_4d_dT =
            Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Zero();
    }

    return ip_constitutive_variables;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::size_t TH2MLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::setIPDataInitialConditions(std::string const& name,
                                                 double const* values,
                                                 int const integration_order)
{
    if (integration_order !=
        static_cast<int>(_integration_method.getIntegrationOrder()))
    {
        OGS_FATAL(
            "Setting integration point initial conditions; The integration "
            "order of the local assembler for element {:d} is different "
            "from the integration order in the initial condition.",
            _element.getID());
    }

    if (name == "sigma_ip")
    {
        if (_process_data.initial_stress != nullptr)
        {
            OGS_FATAL(
                "Setting initial conditions for stress from integration "
                "point data and from a parameter '{:s}' is not possible "
                "simultaneously.",
                _process_data.initial_stress->name);
        }
        return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
            values, _ip_data, &IpData::sigma_eff);
    }

    if (name == "saturation_ip")
    {
        return ProcessLib::setIntegrationPointScalarData(values, _ip_data,
                                                         &IpData::s_L);
    }
    if (name == "epsilon_ip")
    {
        return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
            values, _ip_data, &IpData::eps);
    }
    return 0;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void TH2MLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                        IntegrationMethod, DisplacementDim>::
    setInitialConditionsConcrete(std::vector<double> const& local_x,
                                 double const t,
                                 bool const /*use_monolithic_scheme*/,
                                 int const /*process_id*/)
{
    [[maybe_unused]] auto const matrix_size =
        gas_pressure_size + capillary_pressure_size + temperature_size +
        displacement_size;

    assert(local_x.size() == matrix_size);

    updateConstitutiveVariables(
        Eigen::Map<Eigen::VectorXd const>(local_x.data(), local_x.size()),
        Eigen::VectorXd::Zero(matrix_size), t, 0);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto& ip_data = _ip_data[ip];
        ip_data.pushBackState();
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void TH2MLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::assemble(double const t, double const dt,
                               std::vector<double> const& local_x,
                               std::vector<double> const& local_x_dot,
                               std::vector<double>& local_M_data,
                               std::vector<double>& local_K_data,
                               std::vector<double>& local_rhs_data)
{
    auto const matrix_size = gas_pressure_size + capillary_pressure_size +
                             temperature_size + displacement_size;
    assert(local_x.size() == matrix_size);

    auto const gas_pressure = Eigen::Map<VectorType<gas_pressure_size> const>(
        local_x.data() + gas_pressure_index, gas_pressure_size);

    auto const capillary_pressure =
        Eigen::Map<VectorType<capillary_pressure_size> const>(
            local_x.data() + capillary_pressure_index, capillary_pressure_size);

    auto const capillary_pressure_dot =
        Eigen::Map<VectorType<capillary_pressure_size> const>(
            local_x_dot.data() + capillary_pressure_index,
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

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    updateConstitutiveVariables(
        Eigen::Map<Eigen::VectorXd const>(local_x.data(), local_x.size()),
        Eigen::Map<Eigen::VectorXd const>(local_x_dot.data(),
                                          local_x_dot.size()),
        t, dt);

    for (unsigned int_point = 0; int_point < n_integration_points; int_point++)
    {
        pos.setIntegrationPoint(int_point);
        auto& ip = _ip_data[int_point];

        auto const& Np = ip.N_p;
        auto const& NT = Np;
        auto const& Nu = ip.N_u;

        auto const& NpT = Np.transpose().eval();
        auto const& NTT = NT.transpose().eval();

        auto const& gradNp = ip.dNdx_p;
        auto const& gradNT = gradNp;
        auto const& gradNu = ip.dNdx_u;

        auto const& gradNpT = gradNp.transpose().eval();
        auto const& gradNTT = gradNT.transpose().eval();

        auto const& Nu_op = ip.N_u_op;
        auto const& w = ip.integration_weight;

        auto const& m = Invariants::identity2;

        auto const mT = m.transpose().eval();

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                _element, Nu);

        auto const Bu =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                gradNu, Nu, x_coord, _is_axially_symmetric);

        auto const BuT = Bu.transpose().eval();

        double const pCap = Np.dot(capillary_pressure);

        GlobalDimVectorType const gradpGR = gradNp * gas_pressure;
        GlobalDimVectorType const gradpCap = gradNp * capillary_pressure;

        double const pCap_dot = Np.dot(capillary_pressure_dot);
        auto& beta_T_SR = ip.beta_T_SR;

        auto const I =
            Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Identity();

        const double sD_G = ip.diffusion_coefficient_vapour;
        const double sD_L = ip.diffusion_coefficient_solvate;

        auto const D_C_G = (sD_G * I).eval();
        auto const D_W_G = (sD_G * I).eval();
        auto const D_C_L = (sD_L * I).eval();
        auto const D_W_L = (sD_L * I).eval();

        auto& k_S = ip.k_S;

        auto& s_L = ip.s_L;
        auto const s_G = 1. - s_L;
        auto const s_L_dot = (s_L - ip.s_L_prev) / dt;

        auto& alpha_B = ip.alpha_B;
        auto& beta_p_SR = ip.beta_p_SR;

        auto const& b = _process_data.specific_body_force;

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

        auto const k_over_mu_G = (k_S * ip.k_rel_G / ip.muGR).eval();
        auto const k_over_mu_L = (k_S * ip.k_rel_L / ip.muLR).eval();

        GlobalDimVectorType const w_GS =
            k_over_mu_G * ip.rhoGR * b - k_over_mu_G * gradpGR;

        GlobalDimVectorType const w_LS = k_over_mu_L * gradpCap +
                                         k_over_mu_L * ip.rhoLR * b -
                                         k_over_mu_L * gradpGR;

        // ---------------------------------------------------------------------
        // C-component equation
        // ---------------------------------------------------------------------

        MCpG.noalias() += NpT * rho_C_FR * (alpha_B - phi) * beta_p_SR * Np * w;
        MCpC.noalias() -=
            NpT * rho_C_FR * (alpha_B - phi) * beta_p_SR * s_L * Np * w;

        if (_process_data.apply_mass_lumping)
        {
            if (pCap_dot != 0.)  // avoid division by Zero
            {
                MCpC.noalias() +=
                    NpT *
                    (phi * (ip.rhoCLR - ip.rhoCGR) -
                     rho_C_FR * pCap * (alpha_B - phi) * beta_p_SR) *
                    s_L_dot / pCap_dot * Np * w;
            }
        }

        MCT.noalias() -= NpT * rho_C_FR * (alpha_B - phi) * beta_T_SR * Np * w;
        MCu.noalias() += NpT * rho_C_FR * alpha_B * mT * Bu * w;

        auto const advection_C_G = (ip.rhoCGR * k_over_mu_G).eval();
        auto const advection_C_L = (ip.rhoCLR * k_over_mu_L).eval();
        auto const diffusion_C_G_p =
            (phi_G * ip.rhoGR * D_C_G * ip.dxmCG_dpGR).eval();
        auto const diffusion_C_L_p =
            (phi_L * ip.rhoLR * D_C_L * ip.dxmCL_dpLR).eval();
        auto const diffusion_C_G_T =
            (phi_G * ip.rhoGR * D_C_G * ip.dxmCG_dT).eval();
        auto const diffusion_C_L_T =
            (phi_L * ip.rhoLR * D_C_L * ip.dxmCL_dT).eval();

        auto const advection_C = (advection_C_G + advection_C_L).eval();
        auto const diffusion_C_p = (diffusion_C_G_p + diffusion_C_L_p).eval();
        auto const diffusion_C_T = (diffusion_C_G_T + diffusion_C_L_T).eval();

        LCpG.noalias() += gradNpT * (advection_C + diffusion_C_p) * gradNp * w;

        LCpC.noalias() -=
            gradNpT * (advection_C_L + diffusion_C_L_p) * gradNp * w;

        LCT.noalias() += gradNpT * (diffusion_C_T)*gradNp * w;

        fC.noalias() += gradNpT *
                        (advection_C_G * ip.rhoGR + advection_C_L * ip.rhoLR) *
                        b * w;

        if (!_process_data.apply_mass_lumping)
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

        if (_process_data.apply_mass_lumping)
        {
            if (pCap_dot != 0.)  // avoid division by Zero
            {
                MWpC.noalias() +=
                    NpT *
                    (phi * (ip.rhoWLR - ip.rhoWGR) -
                     rho_W_FR * pCap * (alpha_B - phi) * beta_p_SR) *
                    s_L_dot / pCap_dot * Np * w;
            }
        }

        MWT.noalias() -= NpT * rho_W_FR * (alpha_B - phi) * beta_T_SR * Np * w;

        MWu.noalias() += NpT * rho_W_FR * alpha_B * mT * Bu * w;

        auto const advection_W_G = (ip.rhoWGR * k_over_mu_G).eval();
        auto const advection_W_L = (ip.rhoWLR * k_over_mu_L).eval();
        auto const diffusion_W_G_p =
            (phi_G * ip.rhoGR * D_W_G * ip.dxmWG_dpGR).eval();
        auto const diffusion_W_L_p =
            (phi_L * ip.rhoLR * D_W_L * ip.dxmWL_dpLR).eval();
        auto const diffusion_W_G_T =
            (phi_G * ip.rhoGR * D_W_G * ip.dxmWG_dT).eval();
        auto const diffusion_W_L_T =
            (phi_L * ip.rhoLR * D_W_L * ip.dxmWL_dT).eval();

        auto const advection_W = (advection_W_G + advection_W_L).eval();
        auto const diffusion_W_p = (diffusion_W_G_p + diffusion_W_L_p).eval();
        auto const diffusion_W_T = (diffusion_W_G_T + diffusion_W_L_T).eval();

        LWpG.noalias() += gradNpT * (advection_W + diffusion_W_p) * gradNp * w;

        LWpC.noalias() -=
            gradNpT * (advection_W_L + diffusion_W_L_p) * gradNp * w;

        LWT.noalias() += gradNpT * (diffusion_W_T)*gradNp * w;

        fW.noalias() += gradNpT *
                        (advection_W_G * ip.rhoGR + advection_W_L * ip.rhoLR) *
                        b * w;

        if (!_process_data.apply_mass_lumping)
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
            gradNTT * (ip.rhoGR * h_G * w_GS + ip.rhoLR * h_L * w_LS) * w;

        fT.noalias() +=
            NTT * (ip.rhoGR * w_GS.transpose() + ip.rhoLR * w_LS.transpose()) *
            b * w;

        // ---------------------------------------------------------------------
        //  - displacement equation
        // ---------------------------------------------------------------------

        KUpG.noalias() -= (BuT * alpha_B * m * Np) * w;

        KUpC.noalias() += (BuT * alpha_B * s_L * m * Np) * w;

        fU.noalias() -= (BuT * ip.sigma_eff - Nu_op.transpose() * rho * b) * w;

        if (_process_data.apply_mass_lumping)
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
          typename IntegrationMethod, int DisplacementDim>
void TH2MLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                        IntegrationMethod, DisplacementDim>::
    assembleWithJacobian(double const t, double const dt,
                         std::vector<double> const& local_x,
                         std::vector<double> const& local_xdot,
                         const double /*dxdot_dx*/, const double /*dx_dx*/,
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

    auto const gas_pressure_dot =
        Eigen::Map<VectorType<gas_pressure_size> const>(
            local_xdot.data() + gas_pressure_index, gas_pressure_size);

    auto const capillary_pressure_dot =
        Eigen::Map<VectorType<capillary_pressure_size> const>(
            local_xdot.data() + capillary_pressure_index,
            capillary_pressure_size);

    auto const temperature_dot = Eigen::Map<VectorType<temperature_size> const>(
        local_xdot.data() + temperature_index, temperature_size);

    auto const displacement_dot =
        Eigen::Map<VectorType<displacement_size> const>(
            local_xdot.data() + displacement_index, displacement_size);

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

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    auto const ip_constitutive_variables = updateConstitutiveVariables(
        Eigen::Map<Eigen::VectorXd const>(local_x.data(), local_x.size()),
        Eigen::Map<Eigen::VectorXd const>(local_xdot.data(), local_xdot.size()),
        t, dt);

    for (unsigned int_point = 0; int_point < n_integration_points; int_point++)
    {
        pos.setIntegrationPoint(int_point);
        auto& ip = _ip_data[int_point];
        auto& ip_cv = ip_constitutive_variables[int_point];

        auto const& Np = ip.N_p;
        auto const& NT = Np;
        auto const& Nu = ip.N_u;

        auto const& NpT = Np.transpose().eval();
        auto const& NTT = NT.transpose().eval();

        auto const& gradNp = ip.dNdx_p;
        auto const& gradNT = gradNp;
        auto const& gradNu = ip.dNdx_u;

        auto const& gradNpT = gradNp.transpose().eval();
        auto const& gradNTT = gradNT.transpose().eval();

        auto const& Nu_op = ip.N_u_op;
        auto const& w = ip.integration_weight;

        auto const& m = Invariants::identity2;
        auto const mT = m.transpose().eval();

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                _element, Nu);

        auto const Bu =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                gradNu, Nu, x_coord, _is_axially_symmetric);

        auto const BuT = Bu.transpose().eval();

        double const div_u_dot = Invariants::trace(Bu * displacement_dot);

        double const pCap = Np.dot(capillary_pressure);

        GlobalDimVectorType const gradpGR = gradNp * gas_pressure;
        GlobalDimVectorType const gradpCap = gradNp * capillary_pressure;
        GlobalDimVectorType const gradT = gradNT * temperature;

        double const pGR_dot = Np.dot(gas_pressure_dot);
        double const pCap_dot = Np.dot(capillary_pressure_dot);
        double const T_dot = NT.dot(temperature_dot);
        auto& beta_T_SR = ip.beta_T_SR;

        auto const I =
            Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Identity();

        const double sD_G = ip.diffusion_coefficient_vapour;
        const double sD_L = ip.diffusion_coefficient_solvate;

        auto const D_C_G = (sD_G * I).eval();
        auto const D_W_G = (sD_G * I).eval();
        auto const D_C_L = (sD_L * I).eval();
        auto const D_W_L = (sD_L * I).eval();

        auto& k_S = ip.k_S;

        auto& s_L = ip.s_L;
        auto const s_G = 1. - s_L;
        auto const s_L_dot = (s_L - ip.s_L_prev) / dt;

        auto& alpha_B = ip.alpha_B;
        auto& beta_p_SR = ip.beta_p_SR;

        auto const& b = _process_data.specific_body_force;

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

        auto const k_over_mu_G = (k_S * ip.k_rel_G / ip.muGR).eval();
        auto const k_over_mu_L = (k_S * ip.k_rel_L / ip.muLR).eval();

        GlobalDimVectorType const w_GS =
            k_over_mu_G * ip.rhoGR * b - k_over_mu_G * gradpGR;

        GlobalDimVectorType const w_LS = k_over_mu_L * gradpCap +
                                         k_over_mu_L * ip.rhoLR * b -
                                         k_over_mu_L * gradpGR;

        // ---------------------------------------------------------------------
        // C-component equation
        // ---------------------------------------------------------------------

        MCpG.noalias() += NpT * rho_C_FR * (alpha_B - phi) * beta_p_SR * Np * w;
        MCpC.noalias() -=
            NpT * rho_C_FR * (alpha_B - phi) * beta_p_SR * s_L * Np * w;

        if (_process_data.apply_mass_lumping)
        {
            if (pCap_dot != 0.)  // avoid division by Zero
            {
                MCpC.noalias() +=
                    NpT *
                    (phi * (ip.rhoCLR - ip.rhoCGR) -
                     rho_C_FR * pCap * (alpha_B - phi) * beta_p_SR) *
                    s_L_dot / pCap_dot * Np * w;
            }
        }

        MCT.noalias() -= NpT * rho_C_FR * (alpha_B - phi) * beta_T_SR * Np * w;
        // d (fC_4_MCT * T_dot)/d T
        local_Jac
            .template block<C_size, temperature_size>(C_index,
                                                      temperature_index)
            .noalias() += NpT * ip_cv.dfC_4_MCT_dT * T_dot * NT * w;

        MCu.noalias() += NpT * rho_C_FR * alpha_B * mT * Bu * w;
        // d (fC_4_MCu * u_dot)/d T
        local_Jac
            .template block<C_size, temperature_size>(C_index,
                                                      temperature_index)
            .noalias() += NpT * ip_cv.dfC_4_MCu_dT * div_u_dot * NT * w;

        auto const advection_C_G = (ip.rhoCGR * k_over_mu_G).eval();
        auto const advection_C_L = (ip.rhoCLR * k_over_mu_L).eval();
        auto const diffusion_C_G_p =
            (phi_G * ip.rhoGR * D_C_G * ip.dxmCG_dpGR).eval();
        auto const diffusion_C_L_p =
            (phi_L * ip.rhoLR * D_C_L * ip.dxmCL_dpLR).eval();
        auto const diffusion_C_G_T =
            (phi_G * ip.rhoGR * D_C_G * ip.dxmCG_dT).eval();
        auto const diffusion_C_L_T =
            (phi_L * ip.rhoLR * D_C_L * ip.dxmCL_dT).eval();

        auto const advection_C = (advection_C_G + advection_C_L).eval();
        auto const diffusion_C_p = (diffusion_C_G_p + diffusion_C_L_p).eval();
        auto const diffusion_C_T = (diffusion_C_G_T + diffusion_C_L_T).eval();

        LCpG.noalias() += gradNpT * (advection_C + diffusion_C_p) * gradNp * w;

        // d (fC_4_LCpG * grad p_GR)/d p_GR
        local_Jac.template block<C_size, C_size>(C_index, C_index).noalias() +=
            gradNpT *
            (ip_cv.dadvection_C_dp_GR
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
            NpT * ip_cv.dfC_4_MCpG_dp_GR * pGR_dot * Np * w;

        // d (fC_4_MCpG * p_GR_dot)/d T
        local_Jac
            .template block<C_size, temperature_size>(C_index,
                                                      temperature_index)
            .noalias() += NpT * ip_cv.dfC_4_MCpG_dT * pGR_dot * NT * w;

        LCpC.noalias() -=
            gradNpT * (advection_C_L + diffusion_C_L_p) * gradNp * w;

        LCT.noalias() += gradNpT * diffusion_C_T * gradNp * w;

        // fC_1
        fC.noalias() += gradNpT *
                        (advection_C_G * ip.rhoGR + advection_C_L * ip.rhoLR) *
                        b * w;

        if (!_process_data.apply_mass_lumping)
        {
            // fC_2 = \int a * s_L_dot
            auto const a = phi * (ip.rhoCLR - ip.rhoCGR) -
                           rho_C_FR * pCap * (alpha_B - phi) * beta_p_SR;
            fC.noalias() -= NpT * a * s_L_dot * w;

            // TODO (naumov) Extend for partially saturated media.
            constexpr double ds_L_dp_GR = 0;
            local_Jac.template block<C_size, C_size>(C_index, C_index)
                .noalias() +=
                NpT * (ip_cv.dfC_2a_dp_GR * s_L_dot - a * ds_L_dp_GR / dt) *
                Np * w;

            local_Jac.template block<C_size, W_size>(C_index, W_index)
                .noalias() +=
                NpT *
                (ip_cv.dfC_2a_dp_cap * s_L_dot + a * ip_cv.ds_L_dp_cap / dt) *
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

        if (_process_data.apply_mass_lumping)
        {
            if (pCap_dot != 0.)  // avoid division by Zero
            {
                MWpC.noalias() +=
                    NpT *
                    (phi * (ip.rhoWLR - ip.rhoWGR) -
                     rho_W_FR * pCap * (alpha_B - phi) * beta_p_SR) *
                    s_L_dot / pCap_dot * Np * w;
            }
        }

        MWT.noalias() -= NpT * rho_W_FR * (alpha_B - phi) * beta_T_SR * Np * w;

        MWu.noalias() += NpT * rho_W_FR * alpha_B * mT * Bu * w;

        auto const advection_W_G = (ip.rhoWGR * k_over_mu_G).eval();
        auto const advection_W_L = (ip.rhoWLR * k_over_mu_L).eval();
        auto const diffusion_W_G_p = phi_G * ip.rhoGR * D_W_G * ip.dxmWG_dpGR;
        auto const diffusion_W_L_p = phi_L * ip.rhoLR * D_W_L * ip.dxmWL_dpLR;
        auto const diffusion_W_G_T = phi_G * ip.rhoGR * D_W_G * ip.dxmWG_dT;
        auto const diffusion_W_L_T = phi_L * ip.rhoLR * D_W_L * ip.dxmWL_dT;

        auto const advection_W = advection_W_G + advection_W_L;
        auto const diffusion_W_p = diffusion_W_G_p + diffusion_W_L_p;
        auto const diffusion_W_T = diffusion_W_G_T + diffusion_W_L_T;

        LWpG.noalias() += gradNpT * (advection_W + diffusion_W_p) * gradNp * w;

        // fW_4 LWpG' parts; LWpG = \int grad (a + d) grad
        local_Jac.template block<W_size, C_size>(W_index, C_index).noalias() -=
            gradNpT * (ip_cv.dfW_4a_dp_GR + ip_cv.dfW_4d_dp_GR) * gradpGR * Np *
            w;

        local_Jac.template block<W_size, W_size>(W_index, W_index).noalias() -=
            gradNpT * (ip_cv.dfW_4a_dp_cap + ip_cv.dfW_4d_dp_cap) * gradpGR *
            Np * w;

        local_Jac
            .template block<W_size, temperature_size>(W_index,
                                                      temperature_index)
            .noalias() -=
            gradNpT * (ip_cv.dfW_4a_dT + ip_cv.dfW_4d_dT) * gradpGR * NT * w;

        LWpC.noalias() -=
            gradNpT * (advection_W_L + diffusion_W_L_p) * gradNp * w;

        LWT.noalias() += gradNpT * (diffusion_W_T)*gradNp * w;

        // fW_1
        fW.noalias() += gradNpT *
                        (advection_W_G * ip.rhoGR + advection_W_L * ip.rhoLR) *
                        b * w;

        // fW_2 = \int (f - g) * s_L_dot
        if (!_process_data.apply_mass_lumping)
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
                 (f - g) * ip_cv.ds_L_dp_cap / dt) *
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
        local_Jac
            .template block<temperature_size, C_size>(temperature_index,
                                                      C_index)
            .noalias() += gradNTT * ip_cv.dlambda_dp_GR * gradT * Np * w;

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
            gradNTT * (ip.rhoGR * h_G * w_GS + ip.rhoLR * h_L * w_LS) * w;

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
            NTT * (ip.rhoGR * w_GS.transpose() + ip.rhoLR * w_LS.transpose()) *
            b * w;

        // ---------------------------------------------------------------------
        //  - displacement equation
        // ---------------------------------------------------------------------

        KUpG.noalias() -= (BuT * alpha_B * m * Np) * w;

        // dfU_2/dp_GR part i.e. part of the d(KUpG*p_GR)/dp_GR derivative is
        // dKUpG/dp_GR + KUpG. The former is zero, the latter is handled below.

        KUpC.noalias() += (BuT * alpha_B * s_L * m * Np) * w;

        // dfU_2/dp_LR part i.e. part of the d(KUpC*p_cap)/dp_LR derivative is
        // dKUpC/dp_LR + KUpC. The latter is handled below, the former here:
        local_Jac
            .template block<displacement_size, W_size>(displacement_index,
                                                       W_index)
            .noalias() += BuT * alpha_B * ip_cv.ds_L_dp_cap * pCap * m * Np * w;

        local_Jac
            .template block<displacement_size, displacement_size>(
                displacement_index, displacement_index)
            .noalias() += BuT * ip_cv.C * Bu * w;

        // fU_1
        fU.noalias() -= (BuT * ip.sigma_eff - Nu_op.transpose() * rho * b) * w;

        // KuT
        local_Jac
            .template block<displacement_size, temperature_size>(
                displacement_index, temperature_index)
            .noalias() -= BuT * (ip_cv.C * ip.alpha_T_SR) * NT * w;

        /* TODO (naumov) Test with gravity needed to check this Jacobian part.
        local_Jac
            .template block<displacement_size, temperature_size>(
                displacement_index, temperature_index)
            .noalias() += Nu_op.transpose() * ip_cv.drho_dT * b * Nu_op * w;
            */

        if (_process_data.apply_mass_lumping)
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
                    LCT * temperature + MCpG * gas_pressure_dot +
                    MCpC * capillary_pressure_dot + MCT * temperature_dot +
                    MCu * displacement_dot;

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
                    LWT * temperature + MWpG * gas_pressure_dot +
                    MWpC * capillary_pressure_dot + MWT * temperature_dot +
                    MWu * displacement_dot;

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
    fT.noalias() -= KTT * temperature + MTu * displacement_dot;

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
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const&
TH2MLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                   IntegrationMethod, DisplacementDim>::
    getIntPtDarcyVelocityGas(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const
{
    auto const num_intpts = _ip_data.size();

    constexpr int process_id = 0;  // monolithic scheme;
    auto const indices =
        NumLib::getIndices(_element.getID(), *dof_table[process_id]);
    assert(!indices.empty());
    auto const local_x = x[process_id]->get(indices);

    cache.clear();
    auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, DisplacementDim, num_intpts);

    auto const pGR =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            gas_pressure_size> const>(local_x.data() + gas_pressure_index,
                                      gas_pressure_size);
    auto const pCap =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            capillary_pressure_size> const>(
            local_x.data() + capillary_pressure_index, capillary_pressure_size);
    auto const T =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            temperature_size> const>(local_x.data() + temperature_index,
                                     temperature_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());

    auto const& medium = *_process_data.media_map->getMedium(_element.getID());
    auto const& gas_phase = medium.phase("Gas");

    MPL::VariableArray vars;

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);

        auto const& N_p = _ip_data[ip].N_p;

        vars[static_cast<int>(MPL::Variable::temperature)] =
            N_p.dot(T);  // N_p = N_T
        vars[static_cast<int>(MPL::Variable::phase_pressure)] = N_p.dot(pGR);
        vars[static_cast<int>(MPL::Variable::capillary_pressure)] =
            N_p.dot(pCap);

        // TODO (naumov) Temporary value not used by current material
        // models. Need extension of secondary variables interface.
        double const dt = std::numeric_limits<double>::quiet_NaN();

        auto const mu_GR = gas_phase.property(MPL::PropertyType::viscosity)
                               .template value<double>(vars, pos, t, dt);

        GlobalDimMatrixType k_S = MPL::formEigenTensor<DisplacementDim>(
            medium.property(MPL::PropertyType::permeability)
                .value(vars, pos, t, dt));

        auto const s_L = medium.property(MPL::PropertyType::saturation)
                             .template value<double>(vars, pos, t, dt);

        vars[static_cast<int>(MPL::Variable::liquid_saturation)] = s_L;

        auto const k_rel =
            medium
                .property(
                    MPL::PropertyType::relative_permeability_nonwetting_phase)
                .template value<double>(vars, pos, t, dt);

        auto const k_over_mu = k_S * k_rel / mu_GR;

        vars[static_cast<int>(MPL::Variable::molar_mass)] = 0.1;
        auto const rho_GR = gas_phase.property(MPL::PropertyType::density)
                                .template value<double>(vars, pos, t, dt);
        auto const& b = _process_data.specific_body_force;

        // Compute the velocity
        auto const& dNdx_p = _ip_data[ip].dNdx_p;
        cache_matrix.col(ip).noalias() =
            -k_over_mu * dNdx_p * pGR + k_over_mu * rho_GR * b;
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const&
TH2MLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                   IntegrationMethod, DisplacementDim>::
    getIntPtDarcyVelocityLiquid(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const
{
    auto const num_intpts = _ip_data.size();

    constexpr int process_id = 0;  // monolithic scheme;
    auto const indices =
        NumLib::getIndices(_element.getID(), *dof_table[process_id]);
    assert(!indices.empty());
    auto const local_x = x[process_id]->get(indices);

    cache.clear();
    auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, DisplacementDim, num_intpts);

    auto const pGR =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            gas_pressure_size> const>(local_x.data() + gas_pressure_index,
                                      gas_pressure_size);
    auto const pCap =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            capillary_pressure_size> const>(
            local_x.data() + capillary_pressure_index, capillary_pressure_size);
    auto const pLR = pGR - pCap;
    auto const T =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            temperature_size> const>(local_x.data() + temperature_index,
                                     temperature_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());

    auto const& medium = *_process_data.media_map->getMedium(_element.getID());
    auto const& liquid_phase = medium.phase("AqueousLiquid");

    MPL::VariableArray vars;

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);

        auto const& N_p = _ip_data[ip].N_p;

        vars[static_cast<int>(MPL::Variable::temperature)] = N_p.dot(T);
        vars[static_cast<int>(MPL::Variable::phase_pressure)] = N_p.dot(pGR);
        vars[static_cast<int>(MPL::Variable::liquid_phase_pressure)] =
            N_p.dot(pLR);
        vars[static_cast<int>(MPL::Variable::capillary_pressure)] =
            N_p.dot(pCap);

        // TODO (naumov) Temporary value not used by current material
        // models. Need extension of secondary variables interface.
        double const dt = std::numeric_limits<double>::quiet_NaN();

        auto const mu_LR = liquid_phase.property(MPL::PropertyType::viscosity)
                               .template value<double>(vars, pos, t, dt);
        GlobalDimMatrixType k_S = MPL::formEigenTensor<DisplacementDim>(
            medium.property(MPL::PropertyType::permeability)
                .value(vars, pos, t, dt));

        auto const s_L = medium.property(MPL::PropertyType::saturation)
                             .template value<double>(vars, pos, t, dt);

        vars[static_cast<int>(MPL::Variable::liquid_saturation)] = s_L;

        auto const k_rel =
            medium.property(MPL::PropertyType::relative_permeability)
                .template value<double>(vars, pos, t, dt);

        auto const k_over_mu = k_S * k_rel / mu_LR;

        vars[static_cast<int>(MPL::Variable::molar_fraction)] = 1.0;

        auto const cCL = [&]()
        {
            if (liquid_phase.hasProperty(MPL::PropertyType::concentration))
            {
                return liquid_phase.property(MPL::PropertyType::concentration)
                    .template value<double>(vars, pos, t, dt);  // in mol*m^(-3)
            }
            return 0.;
        }();

        vars[static_cast<int>(MPL::Variable::concentration)] = cCL;

        auto const rho_LR = liquid_phase.property(MPL::PropertyType::density)
                                .template value<double>(vars, pos, t, dt);
        auto const& b = _process_data.specific_body_force;

        // Compute the velocity
        auto const& dNdx_p = _ip_data[ip].dNdx_p;
        cache_matrix.col(ip).noalias() =
            -k_over_mu * dNdx_p * pLR + k_over_mu * rho_LR * b;
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void TH2MLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                        IntegrationMethod, DisplacementDim>::
    computeSecondaryVariableConcrete(double const t, double const dt,
                                     Eigen::VectorXd const& local_x,
                                     Eigen::VectorXd const& local_x_dot)
{
    auto const gas_pressure =
        local_x.template segment<gas_pressure_size>(gas_pressure_index);
    auto const capillary_pressure =
        local_x.template segment<capillary_pressure_size>(
            capillary_pressure_index);
    auto const liquid_pressure = (gas_pressure - capillary_pressure).eval();

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(_element, _is_axially_symmetric, gas_pressure,
                         *_process_data.gas_pressure_interpolated);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(_element, _is_axially_symmetric, capillary_pressure,
                         *_process_data.capillary_pressure_interpolated);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(_element, _is_axially_symmetric, liquid_pressure,
                         *_process_data.liquid_pressure_interpolated);

    auto const temperature =
        local_x.template segment<temperature_size>(temperature_index);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(_element, _is_axially_symmetric, temperature,
                         *_process_data.temperature_interpolated);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    double saturation_avg = 0;

    updateConstitutiveVariables(local_x, local_x_dot, t, dt);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        saturation_avg += _ip_data[ip].s_L;
    }
    saturation_avg /= n_integration_points;
    (*_process_data.element_saturation)[_element.getID()] = saturation_avg;
}

}  // namespace TH2M
}  // namespace ProcessLib
