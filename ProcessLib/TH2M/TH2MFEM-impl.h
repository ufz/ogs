/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "NumLib/Fem/Interpolation.h"
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
    auto const gas_pressure_prev =
        local_x_prev.template segment<gas_pressure_size>(gas_pressure_index);
    auto const capillary_pressure =
        local_x.template segment<capillary_pressure_size>(
            capillary_pressure_index);
    auto const capillary_pressure_prev =
        local_x_prev.template segment<capillary_pressure_size>(
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
    ConstitutiveRelations::MediaData media_data{medium};

    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();

    std::vector<ConstitutiveRelations::ConstitutiveData<DisplacementDim>>
        ip_constitutive_data(n_integration_points);
    std::vector<ConstitutiveRelations::ConstitutiveTempData<DisplacementDim>>
        ip_constitutive_variables(n_integration_points);

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
            std::nullopt, this->element_.getID(),
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
        double const pG = Np.dot(gas_pressure);
        double const pG_prev = Np.dot(gas_pressure_prev);
        double const pCap = Np.dot(capillary_pressure);
        double const pCap_prev = Np.dot(capillary_pressure_prev);
        ConstitutiveRelations::TemperatureData const T_data{T, T_prev};
        ConstitutiveRelations::GasPressureData const pGR_data{pG, pG_prev};
        ConstitutiveRelations::CapillaryPressureData const pCap_data{pCap,
                                                                     pCap_prev};
        ConstitutiveRelations::ReferenceTemperatureData const T0{
            this->process_data_.reference_temperature(t, pos)[0]};
        ConstitutiveRelations::GasPressureGradientData<DisplacementDim> const
            grad_p_GR{gradNp * gas_pressure};
        ConstitutiveRelations::CapillaryPressureGradientData<
            DisplacementDim> const grad_p_cap{gradNp * capillary_pressure};
        ConstitutiveRelations::TemperatureGradientData<DisplacementDim> const
            grad_T{gradNp * temperature};

        // medium properties
        models.elastic_tangent_stiffness_model.eval({pos, t, dt}, T_data,
                                                    ip_cv.C_el_data);

        models.biot_model.eval({pos, t, dt}, media_data, ip_cv.biot_data);

        auto const Bu =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                gradNu, Nu, x_coord, this->is_axially_symmetric_);

        ip_out.eps_data.eps.noalias() = Bu * displacement;
        models.S_L_model.eval({pos, t, dt}, media_data, pCap_data,
                              current_state.S_L_data);

        models.chi_S_L_model.eval({pos, t, dt}, media_data,
                                  current_state.S_L_data,
                                  current_state.chi_S_L);

        models.chi_S_L_prev_model.eval({pos, t, dt}, media_data,
                                       prev_state.S_L_data, prev_state.chi_S_L);

        // solid phase compressibility
        models.beta_p_SR_model.eval({pos, t, dt}, ip_cv.biot_data,
                                    ip_cv.C_el_data, ip_cv.beta_p_SR);

        // If there is swelling stress rate, compute swelling stress.
        models.swelling_model.eval(
            {pos, t, dt}, media_data, ip_cv.C_el_data, current_state.S_L_data,
            prev_state.S_L_data, prev_state.swelling_data,
            current_state.swelling_data, ip_cv.swelling_data);

        // solid phase linear thermal expansion coefficient
        models.s_therm_exp_model.eval({pos, t, dt}, media_data, T_data, T0,
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

        models.total_stress_model.eval(current_state.eff_stress_data,
                                       ip_cv.biot_data, current_state.chi_S_L,
                                       pGR_data, pCap_data,
                                       ip_cv.total_stress_data);

        models.permeability_model.eval(
            {pos, t, dt}, media_data, current_state.S_L_data, pCap_data, T_data,
            ip_cv.total_stress_data, ip_out.eps_data,
            ip_cv.equivalent_plastic_strain_data, ip_out.permeability_data);

        models.pure_liquid_density_model.eval({pos, t, dt}, media_data,
                                              pGR_data, pCap_data, T_data,
                                              current_state.rho_W_LR);

        models.phase_transition_model.eval(
            {pos, t, dt}, media_data, pGR_data, pCap_data, T_data,
            current_state.rho_W_LR, ip_out.fluid_enthalpy_data,
            ip_out.mass_mole_fractions_data, ip_out.fluid_density_data,
            ip_out.vapour_pressure_data, current_state.constituent_density_data,
            ip_cv.phase_transition_data);

        models.viscosity_model.eval({pos, t, dt}, media_data, T_data,
                                    ip_out.mass_mole_fractions_data,
                                    ip_cv.viscosity_data);

        models.porosity_model.eval(
            {pos, t, dt}, media_data, current_state.S_L_data,
            prev_state.S_L_data, pCap_data, pGR_data, current_state.chi_S_L,
            prev_state.chi_S_L, ip_cv.beta_p_SR, ip_out.eps_data,
            Bu * displacement_prev, prev_state.porosity_data,
            current_state.porosity_data);

        if (medium.hasProperty(MPL::PropertyType::transport_porosity))
        {
            models.transport_porosity_model.eval(
                {pos, t, dt}, media_data, current_state.S_L_data,
                prev_state.S_L_data, pCap_data, pGR_data, current_state.chi_S_L,
                prev_state.chi_S_L, ip_cv.beta_p_SR,
                current_state.mechanical_strain_data,
                prev_state.mechanical_strain_data,
                prev_state.transport_porosity_data, current_state.porosity_data,
                current_state.transport_porosity_data);
        }
        else
        {
            current_state.transport_porosity_data.phi =
                current_state.porosity_data.phi;
        }

        models.solid_density_model.eval(
            {pos, t, dt}, media_data, T_data, current_state.eff_stress_data,
            pCap_data, pGR_data, current_state.chi_S_L,
            current_state.porosity_data, ip_out.solid_density_data);

        models.solid_heat_capacity_model.eval({pos, t, dt}, media_data, T_data,
                                              ip_cv.solid_heat_capacity_data);

        models.thermal_conductivity_model.eval(
            {pos, t, dt}, media_data, T_data, current_state.porosity_data,
            current_state.S_L_data, ip_cv.thermal_conductivity_data);

        models.advection_model.eval(current_state.constituent_density_data,
                                    ip_out.permeability_data,
                                    current_state.rho_W_LR,
                                    ip_cv.viscosity_data,
                                    ip_cv.advection_data);

        models.gravity_model.eval(
            ip_out.fluid_density_data,
            current_state.porosity_data,
            current_state.S_L_data,
            ip_out.solid_density_data,
            ConstitutiveRelations::SpecificBodyForceData<DisplacementDim>{
                this->process_data_.specific_body_force},
            ip_cv.volumetric_body_force);

        models.diffusion_velocity_model.eval(grad_p_cap,
                                             grad_p_GR,
                                             ip_out.mass_mole_fractions_data,
                                             ip_cv.phase_transition_data,
                                             current_state.porosity_data,
                                             current_state.S_L_data,
                                             grad_T,
                                             ip_out.diffusion_velocity_data);

        models.solid_enthalpy_model.eval(ip_cv.solid_heat_capacity_data, T_data,
                                         ip_out.solid_enthalpy_data);

        models.internal_energy_model.eval(ip_out.fluid_density_data,
                                          ip_cv.phase_transition_data,
                                          current_state.porosity_data,
                                          current_state.S_L_data,
                                          ip_out.solid_density_data,
                                          ip_out.solid_enthalpy_data,
                                          current_state.internal_energy_data);

        models.effective_volumetric_enthalpy_model.eval(
            ip_out.fluid_density_data,
            ip_out.fluid_enthalpy_data,
            current_state.porosity_data,
            current_state.S_L_data,
            ip_out.solid_density_data,
            ip_out.solid_enthalpy_data,
            ip_cv.effective_volumetric_enthalpy_data);

        models.fC_1_model.eval(ip_cv.advection_data, ip_out.fluid_density_data,
                               ip_cv.fC_1);

        if (!this->process_data_.apply_mass_lumping)
        {
            models.fC_2a_model.eval(ip_cv.biot_data,
                                    pCap_data,
                                    current_state.constituent_density_data,
                                    current_state.porosity_data,
                                    current_state.S_L_data,
                                    ip_cv.beta_p_SR,
                                    ip_cv.fC_2a);
        }
        models.fC_3a_model.eval(dt,
                                current_state.constituent_density_data,
                                prev_state.constituent_density_data,
                                current_state.S_L_data,
                                ip_cv.fC_3a);

        models.fC_4_LCpG_model.eval(ip_cv.advection_data,
                                    ip_out.fluid_density_data,
                                    ip_cv.phase_transition_data,
                                    current_state.porosity_data,
                                    current_state.S_L_data,
                                    ip_cv.fC_4_LCpG);

        models.fC_4_LCpC_model.eval(ip_cv.advection_data,
                                    ip_out.fluid_density_data,
                                    ip_cv.phase_transition_data,
                                    current_state.porosity_data,
                                    current_state.S_L_data,
                                    ip_cv.fC_4_LCpC);

        models.fC_4_LCT_model.eval(ip_out.fluid_density_data,
                                   ip_cv.phase_transition_data,
                                   current_state.porosity_data,
                                   current_state.S_L_data,
                                   ip_cv.fC_4_LCT);

        models.fC_4_MCpG_model.eval(ip_cv.biot_data,
                                    current_state.constituent_density_data,
                                    current_state.porosity_data,
                                    current_state.S_L_data,
                                    ip_cv.beta_p_SR,
                                    ip_cv.fC_4_MCpG);

        models.fC_4_MCpC_model.eval(ip_cv.biot_data,
                                    pCap_data,
                                    current_state.constituent_density_data,
                                    current_state.porosity_data,
                                    prev_state.S_L_data,
                                    current_state.S_L_data,
                                    ip_cv.beta_p_SR,
                                    ip_cv.fC_4_MCpC);

        models.fC_4_MCT_model.eval(ip_cv.biot_data,
                                   current_state.constituent_density_data,
                                   current_state.porosity_data,
                                   current_state.S_L_data,
                                   ip_cv.s_therm_exp_data,
                                   ip_cv.fC_4_MCT);

        models.fC_4_MCu_model.eval(ip_cv.biot_data,
                                   current_state.constituent_density_data,
                                   current_state.S_L_data,
                                   ip_cv.fC_4_MCu);

        models.fW_1_model.eval(ip_cv.advection_data, ip_out.fluid_density_data,
                               ip_cv.fW_1);

        if (!this->process_data_.apply_mass_lumping)
        {
            models.fW_2_model.eval(ip_cv.biot_data,
                                   pCap_data,
                                   current_state.constituent_density_data,
                                   current_state.porosity_data,
                                   current_state.rho_W_LR,
                                   current_state.S_L_data,
                                   ip_cv.beta_p_SR,
                                   ip_cv.fW_2);
        }
        models.fW_3a_model.eval(dt,
                                current_state.constituent_density_data,
                                prev_state.constituent_density_data,
                                prev_state.rho_W_LR,
                                current_state.rho_W_LR,
                                current_state.S_L_data,
                                ip_cv.fW_3a);

        models.fW_4_LWpG_model.eval(ip_cv.advection_data,
                                    ip_out.fluid_density_data,
                                    ip_cv.phase_transition_data,
                                    current_state.porosity_data,
                                    current_state.S_L_data,
                                    ip_cv.fW_4_LWpG);

        models.fW_4_LWpC_model.eval(ip_cv.advection_data,
                                    ip_out.fluid_density_data,
                                    ip_cv.phase_transition_data,
                                    current_state.porosity_data,
                                    current_state.S_L_data,
                                    ip_cv.fW_4_LWpC);

        models.fW_4_LWT_model.eval(ip_out.fluid_density_data,
                                   ip_cv.phase_transition_data,
                                   current_state.porosity_data,
                                   current_state.S_L_data,
                                   ip_cv.fW_4_LWT);

        models.fW_4_MWpG_model.eval(ip_cv.biot_data,
                                    current_state.constituent_density_data,
                                    current_state.porosity_data,
                                    current_state.rho_W_LR,
                                    current_state.S_L_data,
                                    ip_cv.beta_p_SR,
                                    ip_cv.fW_4_MWpG);

        models.fW_4_MWpC_model.eval(ip_cv.biot_data,
                                    pCap_data,
                                    current_state.constituent_density_data,
                                    current_state.porosity_data,
                                    prev_state.S_L_data,
                                    current_state.rho_W_LR,
                                    current_state.S_L_data,
                                    ip_cv.beta_p_SR,
                                    ip_cv.fW_4_MWpC);

        models.fW_4_MWT_model.eval(ip_cv.biot_data,
                                   current_state.constituent_density_data,
                                   current_state.porosity_data,
                                   current_state.rho_W_LR,
                                   current_state.S_L_data,
                                   ip_cv.s_therm_exp_data,
                                   ip_cv.fW_4_MWT);

        models.fW_4_MWu_model.eval(ip_cv.biot_data,
                                   current_state.constituent_density_data,
                                   current_state.rho_W_LR,
                                   current_state.S_L_data,
                                   ip_cv.fW_4_MWu);

        models.fT_1_model.eval(dt,
                               current_state.internal_energy_data,
                               prev_state.internal_energy_data,
                               ip_cv.fT_1);

        // ---------------------------------------------------------------------
        // Derivatives for Jacobian
        // ---------------------------------------------------------------------

        models.darcy_velocity_model.eval(
            grad_p_cap,
            ip_out.fluid_density_data,
            grad_p_GR,
            ip_out.permeability_data,
            ConstitutiveRelations::SpecificBodyForceData<DisplacementDim>{
                this->process_data_.specific_body_force},
            ip_cv.viscosity_data,
            ip_out.darcy_velocity_data);

        models.fT_2_model.eval(ip_out.darcy_velocity_data,
                               ip_out.fluid_density_data,
                               ip_out.fluid_enthalpy_data,
                               ip_cv.fT_2);

        models.fT_3_model.eval(
            current_state.constituent_density_data,
            ip_out.darcy_velocity_data,
            ip_out.diffusion_velocity_data,
            ip_out.fluid_density_data,
            ip_cv.phase_transition_data,
            ConstitutiveRelations::SpecificBodyForceData<DisplacementDim>{
                this->process_data_.specific_body_force},
            ip_cv.fT_3);

        models.fu_2_KupC_model.eval(ip_cv.biot_data, current_state.chi_S_L,
                                    ip_cv.fu_2_KupC);
    }

    return {ip_constitutive_data, ip_constitutive_variables};
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<ConstitutiveRelations::DerivativesData<DisplacementDim>>
TH2MLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                   DisplacementDim>::
    updateConstitutiveVariablesDerivatives(
        Eigen::VectorXd const& local_x, Eigen::VectorXd const& local_x_prev,
        double const t, double const dt,
        std::vector<
            ConstitutiveRelations::ConstitutiveData<DisplacementDim>> const&
            ip_constitutive_data,
        std::vector<
            ConstitutiveRelations::ConstitutiveTempData<DisplacementDim>> const&
            ip_constitutive_variables,
        ConstitutiveRelations::ConstitutiveModels<DisplacementDim> const&
            models)
{
    [[maybe_unused]] auto const matrix_size =
        gas_pressure_size + capillary_pressure_size + temperature_size +
        displacement_size;

    assert(local_x.size() == matrix_size);

    auto const gas_pressure =
        local_x.template segment<gas_pressure_size>(gas_pressure_index);
    auto const gas_pressure_prev =
        local_x_prev.template segment<gas_pressure_size>(gas_pressure_index);
    auto const temperature =
        local_x.template segment<temperature_size>(temperature_index);
    auto const temperature_prev =
        local_x_prev.template segment<temperature_size>(temperature_index);
    auto const displacement_prev =
        local_x_prev.template segment<displacement_size>(displacement_index);

    auto const capillary_pressure =
        local_x.template segment<capillary_pressure_size>(
            capillary_pressure_index);
    auto const capillary_pressure_prev =
        local_x_prev.template segment<capillary_pressure_size>(
            capillary_pressure_index);

    auto const& medium =
        *this->process_data_.media_map.getMedium(this->element_.getID());
    ConstitutiveRelations::MediaData media_data{medium};

    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();

    std::vector<ConstitutiveRelations::DerivativesData<DisplacementDim>>
        ip_d_data(n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& ip_data = _ip_data[ip];
        auto& ip_dd = ip_d_data[ip];
        auto const& ip_cd = ip_constitutive_data[ip];
        auto const& ip_cv = ip_constitutive_variables[ip];
        auto const& ip_out = this->output_data_[ip];
        auto const& current_state = this->current_states_[ip];
        auto const& prev_state = this->prev_states_[ip];

        auto const& Nu = ip_data.N_u;
        auto const& Np = ip_data.N_p;
        auto const& NT = Np;
        auto const& gradNu = ip_data.dNdx_u;
        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                this->element_, Nu);

        ParameterLib::SpatialPosition const pos{
            std::nullopt, this->element_.getID(),
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    this->element_, Nu))};

        double const T = NT.dot(temperature);
        double const T_prev = NT.dot(temperature_prev);
        double const pG = Np.dot(gas_pressure);
        double const pG_prev = Np.dot(gas_pressure_prev);
        double const pCap = Np.dot(capillary_pressure);
        double const pCap_prev = Np.dot(capillary_pressure_prev);
        ConstitutiveRelations::TemperatureData const T_data{T, T_prev};
        ConstitutiveRelations::GasPressureData const pGR_data{pG, pG_prev};
        ConstitutiveRelations::CapillaryPressureData const pCap_data{pCap,
                                                                     pCap_prev};

        auto const Bu =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                gradNu, Nu, x_coord, this->is_axially_symmetric_);

        models.S_L_model.dEval({pos, t, dt}, media_data, pCap_data,
                               ip_dd.dS_L_dp_cap);

        models.advection_model.dEval(current_state.constituent_density_data,
                                     ip_out.permeability_data,
                                     ip_cv.viscosity_data,
                                     ip_dd.dS_L_dp_cap,
                                     ip_cv.phase_transition_data,
                                     ip_dd.advection_d_data);

        models.porosity_model.dEval(
            {pos, t, dt}, media_data, current_state.S_L_data,
            prev_state.S_L_data, pCap_data, pGR_data, current_state.chi_S_L,
            prev_state.chi_S_L, ip_cv.beta_p_SR, ip_out.eps_data,
            Bu * displacement_prev, prev_state.porosity_data,
            ip_dd.porosity_d_data);

        models.thermal_conductivity_model.dEval(
            {pos, t, dt}, media_data, T_data, current_state.porosity_data,
            ip_dd.porosity_d_data, current_state.S_L_data,
            ip_dd.thermal_conductivity_d_data);

        models.solid_density_model.dEval(
            {pos, t, dt}, media_data, T_data, current_state.eff_stress_data,
            pCap_data, pGR_data, current_state.chi_S_L,
            current_state.porosity_data, ip_dd.solid_density_d_data);

        models.internal_energy_model.dEval(
            ip_out.fluid_density_data,
            ip_cv.phase_transition_data,
            current_state.porosity_data,
            ip_dd.porosity_d_data,
            current_state.S_L_data,
            ip_out.solid_density_data,
            ip_dd.solid_density_d_data,
            ip_out.solid_enthalpy_data,
            ip_cv.solid_heat_capacity_data,
            ip_dd.effective_volumetric_internal_energy_d_data);

        models.effective_volumetric_enthalpy_model.dEval(
            ip_out.fluid_density_data,
            ip_out.fluid_enthalpy_data,
            ip_cv.phase_transition_data,
            current_state.porosity_data,
            ip_dd.porosity_d_data,
            current_state.S_L_data,
            ip_out.solid_density_data,
            ip_dd.solid_density_d_data,
            ip_out.solid_enthalpy_data,
            ip_cv.solid_heat_capacity_data,
            ip_dd.effective_volumetric_enthalpy_d_data);
        if (!this->process_data_.apply_mass_lumping)
        {
            models.fC_2a_model.dEval(ip_cv.biot_data,
                                     pCap_data,
                                     current_state.constituent_density_data,
                                     ip_cv.phase_transition_data,
                                     current_state.porosity_data,
                                     ip_dd.porosity_d_data,
                                     current_state.S_L_data,
                                     ip_dd.dS_L_dp_cap,
                                     ip_cv.beta_p_SR,
                                     ip_dd.dfC_2a);
        }
        models.fC_3a_model.dEval(dt,
                                 current_state.constituent_density_data,
                                 prev_state.constituent_density_data,
                                 ip_cv.phase_transition_data,
                                 current_state.S_L_data,
                                 ip_dd.dS_L_dp_cap,
                                 ip_dd.dfC_3a);

        models.fC_4_LCpG_model.dEval(ip_out.permeability_data,
                                     ip_cv.viscosity_data,
                                     ip_cv.phase_transition_data,
                                     ip_dd.advection_d_data,
                                     ip_dd.dfC_4_LCpG);

        models.fC_4_LCpC_model.dEval(current_state.constituent_density_data,
                                     ip_out.permeability_data,
                                     ip_cv.phase_transition_data,
                                     ip_dd.dS_L_dp_cap,
                                     ip_cv.viscosity_data,
                                     ip_dd.dfC_4_LCpC);

        models.fC_4_MCpG_model.dEval(ip_cv.biot_data,
                                     current_state.constituent_density_data,
                                     ip_cv.phase_transition_data,
                                     current_state.porosity_data,
                                     ip_dd.porosity_d_data,
                                     current_state.S_L_data,
                                     ip_cv.beta_p_SR,
                                     ip_dd.dfC_4_MCpG);

        models.fC_4_MCT_model.dEval(ip_cv.biot_data,
                                    current_state.constituent_density_data,
                                    ip_cv.phase_transition_data,
                                    current_state.porosity_data,
                                    ip_dd.porosity_d_data,
                                    current_state.S_L_data,
                                    ip_cv.s_therm_exp_data,
                                    ip_dd.dfC_4_MCT);

        models.fC_4_MCu_model.dEval(ip_cv.biot_data,
                                    ip_cv.phase_transition_data,
                                    current_state.S_L_data,
                                    ip_dd.dfC_4_MCu);

        if (!this->process_data_.apply_mass_lumping)
        {
            models.fW_2_model.dEval(ip_cv.biot_data,
                                    pCap_data,
                                    current_state.constituent_density_data,
                                    ip_cv.phase_transition_data,
                                    current_state.porosity_data,
                                    ip_dd.porosity_d_data,
                                    current_state.rho_W_LR,
                                    current_state.S_L_data,
                                    ip_dd.dS_L_dp_cap,
                                    ip_cv.beta_p_SR,
                                    ip_dd.dfW_2);
        }

        models.fW_3a_model.dEval(dt,
                                 current_state.constituent_density_data,
                                 ip_cv.phase_transition_data,
                                 prev_state.constituent_density_data,
                                 prev_state.rho_W_LR,
                                 current_state.rho_W_LR,
                                 current_state.S_L_data,
                                 ip_dd.dS_L_dp_cap,
                                 ip_dd.dfW_3a);

        models.fW_4_LWpG_model.dEval(current_state.constituent_density_data,
                                     ip_out.permeability_data,
                                     ip_cv.phase_transition_data,
                                     current_state.rho_W_LR,
                                     ip_dd.dS_L_dp_cap,
                                     ip_cv.viscosity_data,
                                     ip_dd.dfW_4_LWpG);

        models.fW_4_LWpC_model.dEval(ip_cv.advection_data,
                                     ip_out.fluid_density_data,
                                     ip_out.permeability_data,
                                     ip_cv.phase_transition_data,
                                     current_state.porosity_data,
                                     current_state.rho_W_LR,
                                     current_state.S_L_data,
                                     ip_dd.dS_L_dp_cap,
                                     ip_cv.viscosity_data,
                                     ip_dd.dfW_4_LWpC);

        models.fT_1_model.dEval(
            dt, ip_dd.effective_volumetric_internal_energy_d_data, ip_dd.dfT_1);

        models.fT_2_model.dEval(
            ip_out.darcy_velocity_data,
            ip_out.fluid_density_data,
            ip_out.fluid_enthalpy_data,
            ip_out.permeability_data,
            ip_cv.phase_transition_data,
            ConstitutiveRelations::SpecificBodyForceData<DisplacementDim>{
                this->process_data_.specific_body_force},
            ip_cv.viscosity_data,
            ip_dd.dfT_2);

        models.fu_1_KuT_model.dEval(ip_cd.s_mech_data, ip_cv.s_therm_exp_data,
                                    ip_dd.dfu_1_KuT);

        models.fu_2_KupC_model.dEval(ip_cv.biot_data,
                                     current_state.chi_S_L,
                                     pCap_data,
                                     ip_dd.dS_L_dp_cap,
                                     ip_dd.dfu_2_KupC);
    }

    return ip_d_data;
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
        this->solid_material_, *this->process_data_.phase_transition_model_};

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
            std::nullopt, this->element_.getID(),
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    this->element_, ip_data.N_u))};

        double const pCap = Np.dot(capillary_pressure);
        vars.capillary_pressure = pCap;

        double const T = NT.dot(temperature);
        ConstitutiveRelations::TemperatureData const T_data{
            T, T};  // T_prev = T in initialization.
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

            this->current_states_[ip].eff_stress_data.sigma_eff.noalias() +=
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
        this->material_states_[ip].pushBackState();
        this->prev_states_[ip] = this->current_states_[ip];
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

    auto KUU = local_K.template block<displacement_size, displacement_size>(
        displacement_index, displacement_index);

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
        this->solid_material_, *this->process_data_.phase_transition_model_};

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
        auto& ip_cd = ip_constitutive_data[int_point];

        auto& current_state = this->current_states_[int_point];
        auto const& prev_state = this->prev_states_[int_point];

        auto const& Np = ip.N_p;
        auto const& NT = Np;
        auto const& Nu = ip.N_u;
        ParameterLib::SpatialPosition const pos{
            std::nullopt, this->element_.getID(),
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

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                this->element_, Nu);

        auto const Bu =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                gradNu, Nu, x_coord, this->is_axially_symmetric_);

        auto const NTN = (Np.transpose() * Np).eval();
        auto const BTI2N = (Bu.transpose() * Invariants::identity2 * Np).eval();

        double const pCap = Np.dot(capillary_pressure);
        double const pCap_prev = Np.dot(capillary_pressure_prev);

        auto const s_L = current_state.S_L_data.S_L;
        auto const s_L_dot = (s_L - prev_state.S_L_data->S_L) / dt;

        auto const& b = this->process_data_.specific_body_force;

        // ---------------------------------------------------------------------
        // C-component equation
        // ---------------------------------------------------------------------

        MCpG.noalias() += NTN * (ip_cv.fC_4_MCpG.m * w);
        MCpC.noalias() += NTN * (ip_cv.fC_4_MCpC.m * w);

        if (this->process_data_.apply_mass_lumping)
        {
            if (pCap - pCap_prev != 0.)  // avoid division by Zero
            {
                MCpC.noalias() +=
                    NTN * (ip_cv.fC_4_MCpC.ml / (pCap - pCap_prev) * w);
            }
        }

        MCT.noalias() += NTN * (ip_cv.fC_4_MCT.m * w);
        MCu.noalias() += BTI2N.transpose() * (ip_cv.fC_4_MCu.m * w);

        LCpG.noalias() += gradNpT * ip_cv.fC_4_LCpG.L * gradNp * w;

        LCpC.noalias() += gradNpT * ip_cv.fC_4_LCpC.L * gradNp * w;

        LCT.noalias() += gradNpT * ip_cv.fC_4_LCT.L * gradNp * w;

        fC.noalias() += gradNpT * ip_cv.fC_1.A * b * w;

        if (!this->process_data_.apply_mass_lumping)
        {
            fC.noalias() -= NpT * (ip_cv.fC_2a.a * s_L_dot * w);
        }
        // fC_III
        fC.noalias() -=
            NpT * (current_state.porosity_data.phi * ip_cv.fC_3a.a * w);

        // ---------------------------------------------------------------------
        // W-component equation
        // ---------------------------------------------------------------------

        MWpG.noalias() += NTN * (ip_cv.fW_4_MWpG.m * w);
        MWpC.noalias() += NTN * (ip_cv.fW_4_MWpC.m * w);

        if (this->process_data_.apply_mass_lumping)
        {
            if (pCap - pCap_prev != 0.)  // avoid division by Zero
            {
                MWpC.noalias() +=
                    NTN * (ip_cv.fW_4_MWpC.ml / (pCap - pCap_prev) * w);
            }
        }

        MWT.noalias() += NTN * (ip_cv.fW_4_MWT.m * w);

        MWu.noalias() += BTI2N.transpose() * (ip_cv.fW_4_MWu.m * w);

        LWpG.noalias() += gradNpT * ip_cv.fW_4_LWpG.L * gradNp * w;

        LWpC.noalias() += gradNpT * ip_cv.fW_4_LWpC.L * gradNp * w;

        LWT.noalias() += gradNpT * ip_cv.fW_4_LWT.L * gradNp * w;

        fW.noalias() += gradNpT * ip_cv.fW_1.A * b * w;

        if (!this->process_data_.apply_mass_lumping)
        {
            fW.noalias() -= NpT * (ip_cv.fW_2.a * s_L_dot * w);
        }

        fW.noalias() -=
            NpT * (current_state.porosity_data.phi * ip_cv.fW_3a.a * w);

        // ---------------------------------------------------------------------
        //  - temperature equation
        // ---------------------------------------------------------------------

        MTu.noalias() +=
            BTI2N.transpose() *
            (ip_cv.effective_volumetric_enthalpy_data.rho_h_eff * w);

        KTT.noalias() +=
            gradNTT * ip_cv.thermal_conductivity_data.lambda * gradNT * w;

        fT.noalias() -= NTT * (ip_cv.fT_1.m * w);

        fT.noalias() += gradNTT * ip_cv.fT_2.A * w;

        fT.noalias() += gradNTT * ip_cv.fT_3.gradN * w;

        fT.noalias() += NTT * (ip_cv.fT_3.N * w);

        // ---------------------------------------------------------------------
        //  - displacement equation
        // ---------------------------------------------------------------------

        KUpG.noalias() -= BTI2N * (ip_cv.biot_data() * w);

        KUpC.noalias() += BTI2N * (ip_cv.fu_2_KupC.m * w);

        KUU.noalias() +=
            Bu.transpose() * ip_cd.s_mech_data.stiffness_tensor * Bu * w;

        fU.noalias() -=
            (Bu.transpose() * current_state.eff_stress_data.sigma_eff -
             N_u_op(Nu).transpose() * ip_cv.volumetric_body_force()) *
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
        this->solid_material_, *this->process_data_.phase_transition_model_};

    auto const [ip_constitutive_data, ip_constitutive_variables] =
        updateConstitutiveVariables(
            Eigen::Map<Eigen::VectorXd const>(local_x.data(), local_x.size()),
            Eigen::Map<Eigen::VectorXd const>(local_x_prev.data(),
                                              local_x_prev.size()),
            t, dt, models);

    auto const ip_d_data = updateConstitutiveVariablesDerivatives(
        Eigen::Map<Eigen::VectorXd const>(local_x.data(), local_x.size()),
        Eigen::Map<Eigen::VectorXd const>(local_x_prev.data(),
                                          local_x_prev.size()),
        t, dt, ip_constitutive_data, ip_constitutive_variables, models);

    for (unsigned int_point = 0; int_point < n_integration_points; int_point++)
    {
        auto& ip = _ip_data[int_point];
        auto& ip_cd = ip_constitutive_data[int_point];
        auto& ip_dd = ip_d_data[int_point];
        auto& ip_cv = ip_constitutive_variables[int_point];
        auto& current_state = this->current_states_[int_point];
        auto& prev_state = this->prev_states_[int_point];

        auto const& Np = ip.N_p;
        auto const& NT = Np;
        auto const& Nu = ip.N_u;
        ParameterLib::SpatialPosition const pos{
            std::nullopt, this->element_.getID(),
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

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                this->element_, Nu);

        auto const Bu =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                gradNu, Nu, x_coord, this->is_axially_symmetric_);

        auto const NTN = (Np.transpose() * Np).eval();
        auto const BTI2N = (Bu.transpose() * Invariants::identity2 * Np).eval();

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

        auto const& s_L = current_state.S_L_data.S_L;
        auto const s_L_dot = (s_L - prev_state.S_L_data->S_L) / dt;

        auto const& b = this->process_data_.specific_body_force;

        // ---------------------------------------------------------------------
        // C-component equation
        // ---------------------------------------------------------------------

        MCpG.noalias() += NTN * (ip_cv.fC_4_MCpG.m * w);
        MCpC.noalias() += NTN * (ip_cv.fC_4_MCpC.m * w);

        if (this->process_data_.apply_mass_lumping)
        {
            if (pCap - pCap_prev != 0.)  // avoid division by Zero
            {
                MCpC.noalias() +=
                    NTN * (ip_cv.fC_4_MCpC.ml / (pCap - pCap_prev) * w);
            }
        }

        MCT.noalias() += NTN * (ip_cv.fC_4_MCT.m * w);
        // d (fC_4_MCT * T_dot)/d T
        local_Jac
            .template block<C_size, temperature_size>(C_index,
                                                      temperature_index)
            .noalias() += NTN * (ip_dd.dfC_4_MCT.dT * (T - T_prev) / dt * w);

        MCu.noalias() += BTI2N.transpose() * (ip_cv.fC_4_MCu.m * w);
        // d (fC_4_MCu * u_dot)/d T
        local_Jac
            .template block<C_size, temperature_size>(C_index,
                                                      temperature_index)
            .noalias() += NTN * (ip_dd.dfC_4_MCu.dT * div_u_dot * w);

        LCpG.noalias() += gradNpT * ip_cv.fC_4_LCpG.L * gradNp * w;

        // d (fC_4_LCpG * grad p_GR)/d p_GR
        local_Jac.template block<C_size, C_size>(C_index, C_index).noalias() +=
            gradNpT * ip_dd.dfC_4_LCpG.dp_GR * gradpGR * Np * w;

        // d (fC_4_LCpG * grad p_GR)/d p_cap
        local_Jac.template block<C_size, W_size>(C_index, W_index).noalias() +=
            gradNpT * ip_dd.dfC_4_LCpG.dp_cap * gradpGR * Np * w;

        // d (fC_4_LCpG * grad p_GR)/d T
        local_Jac
            .template block<C_size, temperature_size>(C_index,
                                                      temperature_index)
            .noalias() += gradNpT * ip_dd.dfC_4_LCpG.dT * gradpGR * NT * w;

        // d (fC_4_MCpG * p_GR_dot)/d p_GR
        local_Jac.template block<C_size, C_size>(C_index, C_index).noalias() +=
            NTN * (ip_dd.dfC_4_MCpG.dp_GR * (pGR - pGR_prev) / dt * w);

        // d (fC_4_MCpG * p_GR_dot)/d T
        local_Jac
            .template block<C_size, temperature_size>(C_index,
                                                      temperature_index)
            .noalias() +=
            NTN * (ip_dd.dfC_4_MCpG.dT * (pGR - pGR_prev) / dt * w);

        LCpC.noalias() -= gradNpT * ip_cv.fC_4_LCpC.L * gradNp * w;

        /* TODO (naumov) This part is not tested by any of the current ctests.
        // d (fC_4_LCpC * grad p_cap)/d p_GR
        local_Jac.template block<C_size, C_size>(C_index, C_index).noalias() +=
            gradNpT * ip_dd.dfC_4_LCpC.dp_GR * gradpCap * Np * w;
        // d (fC_4_LCpC * grad p_cap)/d p_cap
        local_Jac.template block<C_size, W_size>(C_index, W_index).noalias() +=
            gradNpT * ip_dd.dfC_4_LCpC.dp_cap * gradpCap * Np * w;

        local_Jac
            .template block<C_size, temperature_size>(C_index,
                                                      temperature_index)
            .noalias() += gradNpT * ip_dd.dfC_4_LCpC.dT * gradpCap * Np * w;
        */

        LCT.noalias() += gradNpT * ip_cv.fC_4_LCT.L * gradNp * w;

        // fC_1
        fC.noalias() += gradNpT * ip_cv.fC_1.A * b * w;

        if (!this->process_data_.apply_mass_lumping)
        {
            // fC_2 = \int a * s_L_dot
            fC.noalias() -= NpT * (ip_cv.fC_2a.a * s_L_dot * w);

            local_Jac.template block<C_size, C_size>(C_index, C_index)
                .noalias() +=
                NTN * ((ip_dd.dfC_2a.dp_GR * s_L_dot
                        /*- ip_cv.fC_2a.a * (ds_L_dp_GR = 0) / dt*/) *
                       w);

            local_Jac.template block<C_size, W_size>(C_index, W_index)
                .noalias() +=
                NTN * ((ip_dd.dfC_2a.dp_cap * s_L_dot +
                        ip_cv.fC_2a.a * ip_dd.dS_L_dp_cap() / dt) *
                       w);

            local_Jac
                .template block<C_size, temperature_size>(C_index,
                                                          temperature_index)
                .noalias() += NTN * (ip_dd.dfC_2a.dT * s_L_dot * w);
        }
        {
            // fC_3 = \int phi * a
            fC.noalias() -=
                NpT * (current_state.porosity_data.phi * ip_cv.fC_3a.a * w);

            local_Jac.template block<C_size, C_size>(C_index, C_index)
                .noalias() += NTN * (current_state.porosity_data.phi *
                                     ip_dd.dfC_3a.dp_GR * w);

            local_Jac.template block<C_size, W_size>(C_index, W_index)
                .noalias() += NTN * (current_state.porosity_data.phi *
                                     ip_dd.dfC_3a.dp_cap * w);

            local_Jac
                .template block<C_size, temperature_size>(C_index,
                                                          temperature_index)
                .noalias() +=
                NTN * ((ip_dd.porosity_d_data.dphi_dT * ip_cv.fC_3a.a +
                        current_state.porosity_data.phi * ip_dd.dfC_3a.dT) *
                       w);
        }
        // ---------------------------------------------------------------------
        // W-component equation
        // ---------------------------------------------------------------------

        MWpG.noalias() += NTN * (ip_cv.fW_4_MWpG.m * w);
        MWpC.noalias() += NTN * (ip_cv.fW_4_MWpC.m * w);

        if (this->process_data_.apply_mass_lumping)
        {
            if (pCap - pCap_prev != 0.)  // avoid division by Zero
            {
                MWpC.noalias() +=
                    NTN * (ip_cv.fW_4_MWpC.ml / (pCap - pCap_prev) * w);
            }
        }

        MWT.noalias() += NTN * (ip_cv.fW_4_MWT.m * w);

        MWu.noalias() += BTI2N.transpose() * (ip_cv.fW_4_MWu.m * w);

        LWpG.noalias() += gradNpT * ip_cv.fW_4_LWpG.L * gradNp * w;

        // fW_4 LWpG' parts; LWpG = \int grad (a + d) grad
        local_Jac.template block<W_size, C_size>(W_index, C_index).noalias() +=
            gradNpT * ip_dd.dfW_4_LWpG.dp_GR * gradpGR * Np * w;

        local_Jac.template block<W_size, W_size>(W_index, W_index).noalias() +=
            gradNpT * ip_dd.dfW_4_LWpG.dp_cap * gradpGR * Np * w;

        local_Jac
            .template block<W_size, temperature_size>(W_index,
                                                      temperature_index)
            .noalias() += gradNpT * ip_dd.dfW_4_LWpG.dT * gradpGR * NT * w;

        LWpC.noalias() += gradNpT * ip_cv.fW_4_LWpC.L * gradNp * w;

        // fW_4 LWp_cap' parts; LWpC = \int grad (a + d) grad
        local_Jac.template block<W_size, C_size>(W_index, C_index).noalias() -=
            gradNpT * ip_dd.dfW_4_LWpC.dp_GR * gradpCap * Np * w;

        local_Jac.template block<W_size, W_size>(W_index, W_index).noalias() -=
            gradNpT * ip_dd.dfW_4_LWpC.dp_cap * gradpCap * Np * w;

        local_Jac
            .template block<W_size, temperature_size>(W_index,
                                                      temperature_index)
            .noalias() -= gradNpT * ip_dd.dfW_4_LWpC.dT * gradpCap * NT * w;

        LWT.noalias() += gradNpT * ip_cv.fW_4_LWT.L * gradNp * w;

        // fW_1
        fW.noalias() += gradNpT * ip_cv.fW_1.A * b * w;

        // fW_2 = \int a * s_L_dot
        if (!this->process_data_.apply_mass_lumping)
        {
            fW.noalias() -= NpT * (ip_cv.fW_2.a * s_L_dot * w);

            local_Jac.template block<W_size, C_size>(W_index, C_index)
                .noalias() += NTN * (ip_dd.dfW_2.dp_GR * s_L_dot * w);

            // sign negated because of dp_cap = -dp_LR
            // TODO (naumov) Had to change the sign to get equal Jacobian WW
            // blocks in A2 Test. Where is the error?
            local_Jac.template block<W_size, W_size>(W_index, W_index)
                .noalias() += NTN * ((ip_dd.dfW_2.dp_cap * s_L_dot +
                                      ip_cv.fW_2.a * ip_dd.dS_L_dp_cap() / dt) *
                                     w);

            local_Jac
                .template block<W_size, temperature_size>(W_index,
                                                          temperature_index)
                .noalias() += NTN * (ip_dd.dfW_2.dT * s_L_dot * w);
        }

        // fW_3 = \int phi * a
        fW.noalias() -=
            NpT * (current_state.porosity_data.phi * ip_cv.fW_3a.a * w);

        local_Jac.template block<W_size, C_size>(W_index, C_index).noalias() +=
            NTN * (current_state.porosity_data.phi * ip_dd.dfW_3a.dp_GR * w);

        local_Jac.template block<W_size, W_size>(W_index, W_index).noalias() +=
            NTN * (current_state.porosity_data.phi * ip_dd.dfW_3a.dp_cap * w);

        local_Jac
            .template block<W_size, temperature_size>(W_index,
                                                      temperature_index)
            .noalias() +=
            NTN * ((ip_dd.porosity_d_data.dphi_dT * ip_cv.fW_3a.a +
                    current_state.porosity_data.phi * ip_dd.dfW_3a.dT) *
                   w);

        // ---------------------------------------------------------------------
        //  - temperature equation
        // ---------------------------------------------------------------------

        MTu.noalias() +=
            BTI2N.transpose() *
            (ip_cv.effective_volumetric_enthalpy_data.rho_h_eff * w);

        // dfT_4/dp_GR
        // d (MTu * u_dot)/dp_GR
        local_Jac
            .template block<temperature_size, C_size>(temperature_index,
                                                      C_index)
            .noalias() +=
            NTN * (ip_dd.effective_volumetric_enthalpy_d_data.drho_h_eff_dp_GR *
                   div_u_dot * w);

        // dfT_4/dp_cap
        // d (MTu * u_dot)/dp_cap
        local_Jac
            .template block<temperature_size, W_size>(temperature_index,
                                                      W_index)
            .noalias() -=
            NTN *
            (ip_dd.effective_volumetric_enthalpy_d_data.drho_h_eff_dp_cap *
             div_u_dot * w);

        // dfT_4/dT
        // d (MTu * u_dot)/dT
        local_Jac
            .template block<temperature_size, temperature_size>(
                temperature_index, temperature_index)
            .noalias() +=
            NTN * (ip_dd.effective_volumetric_enthalpy_d_data.drho_h_eff_dT *
                   div_u_dot * w);

        KTT.noalias() +=
            gradNTT * ip_cv.thermal_conductivity_data.lambda * gradNT * w;

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
            .noalias() += gradNTT *
                          ip_dd.thermal_conductivity_d_data.dlambda_dp_cap *
                          gradT * Np * w;

        // d KTT/dT * T
        local_Jac
            .template block<temperature_size, temperature_size>(
                temperature_index, temperature_index)
            .noalias() += gradNTT *
                          ip_dd.thermal_conductivity_d_data.dlambda_dT * gradT *
                          NT * w;

        // fT_1
        fT.noalias() -= NTT * (ip_cv.fT_1.m * w);

        // dfT_1/dp_GR
        local_Jac
            .template block<temperature_size, C_size>(temperature_index,
                                                      C_index)
            .noalias() += NTN * (ip_dd.dfT_1.dp_GR * w);

        // dfT_1/dp_cap
        local_Jac
            .template block<temperature_size, W_size>(temperature_index,
                                                      W_index)
            .noalias() += NTN * (ip_dd.dfT_1.dp_cap * w);

        // dfT_1/dT
        // MTT
        local_Jac
            .template block<temperature_size, temperature_size>(
                temperature_index, temperature_index)
            .noalias() += NTN * (ip_dd.dfT_1.dT * w);

        // fT_2
        fT.noalias() += gradNTT * ip_cv.fT_2.A * w;

        // dfT_2/dp_GR
        local_Jac
            .template block<temperature_size, C_size>(temperature_index,
                                                      C_index)
            .noalias() -=
            // dfT_2/dp_GR first part
            gradNTT * ip_dd.dfT_2.dp_GR_Npart * Np * w +
            // dfT_2/dp_GR second part
            gradNTT * ip_dd.dfT_2.dp_GR_gradNpart * gradNp * w;

        // dfT_2/dp_cap
        local_Jac
            .template block<temperature_size, W_size>(temperature_index,
                                                      W_index)
            .noalias() -=
            // first part of dfT_2/dp_cap
            gradNTT * (-ip_dd.dfT_2.dp_cap_Npart) * Np * w +
            // second part of dfT_2/dp_cap
            gradNTT * (-ip_dd.dfT_2.dp_cap_gradNpart) * gradNp * w;

        // dfT_2/dT
        local_Jac
            .template block<temperature_size, temperature_size>(
                temperature_index, temperature_index)
            .noalias() -= gradNTT * ip_dd.dfT_2.dT * NT * w;

        // fT_3
        fT.noalias() += NTT * (ip_cv.fT_3.N * w);

        fT.noalias() += gradNTT * ip_cv.fT_3.gradN * w;

        // ---------------------------------------------------------------------
        //  - displacement equation
        // ---------------------------------------------------------------------

        KUpG.noalias() -= BTI2N * (ip_cv.biot_data() * w);

        // dfU_2/dp_GR = dKUpG/dp_GR * p_GR + KUpG. The former is zero, the
        // latter is handled below.

        KUpC.noalias() += BTI2N * (ip_cv.fu_2_KupC.m * w);

        // dfU_2/dp_cap = dKUpC/dp_cap * p_cap + KUpC. The former is handled
        // here, the latter below.
        local_Jac
            .template block<displacement_size, W_size>(displacement_index,
                                                       W_index)
            .noalias() += BTI2N * (ip_dd.dfu_2_KupC.dp_cap * w);

        local_Jac
            .template block<displacement_size, displacement_size>(
                displacement_index, displacement_index)
            .noalias() +=
            Bu.transpose() * ip_cd.s_mech_data.stiffness_tensor * Bu * w;

        // fU_1
        fU.noalias() -=
            (Bu.transpose() * current_state.eff_stress_data.sigma_eff -
             N_u_op(Nu).transpose() * ip_cv.volumetric_body_force()) *
            w;

        // KuT
        local_Jac
            .template block<displacement_size, temperature_size>(
                displacement_index, temperature_index)
            .noalias() -= Bu.transpose() * ip_dd.dfu_1_KuT.dT * NT * w;

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

    ConstitutiveRelations::ConstitutiveModels<DisplacementDim> const models{
        this->solid_material_, *this->process_data_.phase_transition_model_};

    updateConstitutiveVariables(local_x, local_x_prev, t, dt, models);
}

}  // namespace TH2M
}  // namespace ProcessLib
