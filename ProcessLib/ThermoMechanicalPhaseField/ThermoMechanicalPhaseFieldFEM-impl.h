/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   ThermoMechanicalPhaseFieldFEM-impl.h
 *  Created on January 8, 2018, 3:00 PM
 */
#pragma once

#include "ThermoMechanicalPhaseFieldFEM.h"

#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
namespace ThermoMechanicalPhaseField
{
template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void ThermoMechanicalPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                DisplacementDim>::
    assembleWithJacobian(double const t, std::vector<double> const& local_x,
                         std::vector<double> const& local_xdot,
                         const double /*dxdot_dx*/, const double /*dx_dx*/,
                         std::vector<double>& /*local_M_data*/,
                         std::vector<double>& /*local_K_data*/,
                         std::vector<double>& local_rhs_data,
                         std::vector<double>& local_Jac_data)
{
    using RhsVector =
        typename ShapeMatricesType::template VectorType<phasefield_size +
                                                        displacement_size>;
    using JacobianMatrix = typename ShapeMatricesType::template MatrixType<
        phasefield_size + displacement_size,
        phasefield_size + displacement_size>;

    auto const local_matrix_size = local_x.size();
    assert(local_matrix_size ==
           temperature_size + phasefield_size + displacement_size);

    auto T = Eigen::Map<typename ShapeMatricesType::template VectorType<
        temperature_size> const>(local_x.data() + temperature_index,
                                 temperature_size);

    auto d = Eigen::Map<
        typename ShapeMatricesType::template VectorType<phasefield_size> const>(
        local_x.data() + phasefield_index, phasefield_size);

    auto u = Eigen::Map<typename ShapeMatricesType::template VectorType<
        displacement_size> const>(local_x.data() + displacement_index,
                                  displacement_size);

    auto T_dot = Eigen::Map<typename ShapeMatricesType::template VectorType<
        temperature_size> const>(local_xdot.data() + temperature_index,
                                 temperature_size);

    auto d_dot = Eigen::Map<typename ShapeMatricesType::template VectorType<
        temperature_size> const>(local_xdot.data() + phasefield_index,
                                 phasefield_size);

    auto local_Jac = MathLib::createZeroedMatrix<JacobianMatrix>(
        local_Jac_data, local_matrix_size, local_matrix_size);

    auto local_rhs = MathLib::createZeroedVector<RhsVector>(local_rhs_data,
                                                            local_matrix_size);

    typename ShapeMatricesType::template MatrixType<displacement_size,
                                                    phasefield_size>
        Kud;
    Kud.setZero(displacement_size, phasefield_size);

    typename ShapeMatricesType::template MatrixType<displacement_size,
                                                    temperature_size>
        KuT;
    KuT.setZero(displacement_size, temperature_size);

    typename ShapeMatricesType::template MatrixType<phasefield_size,
                                                    displacement_size>
        Kdu;
    Kdu.setZero(phasefield_size, displacement_size);

    typename ShapeMatricesType::template MatrixType<phasefield_size,
                                                    temperature_size>
        KdT;
    KdT.setZero(phasefield_size, temperature_size);

    typename ShapeMatricesType::template MatrixType<temperature_size,
                                                    phasefield_size>
        KTd;
    KTd.setZero(temperature_size, phasefield_size);

    typename ShapeMatricesType::NodalMatrixType KTT;
    KTT.setZero(temperature_size, temperature_size);

    typename ShapeMatricesType::NodalMatrixType DTT;
    DTT.setZero(temperature_size, temperature_size);

    typename ShapeMatricesType::NodalMatrixType Kdd;
    Kdd.setZero(phasefield_size, phasefield_size);

    typename ShapeMatricesType::NodalMatrixType Ddd;
    Ddd.setZero(phasefield_size, phasefield_size);

    double const& dt = _process_data.dt;

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;

        auto const& dNdx = _ip_data[ip].dNdx;
        auto const& N = _ip_data[ip].N;

        auto const x_coord =
            interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(_element,
                                                                     N);
        auto const& B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunction::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx, N, x_coord, _is_axially_symmetric);

        auto& eps_m = _ip_data[ip].eps_m;
        auto& eps = _ip_data[ip].eps;

        auto const& C_tensile = _ip_data[ip].C_tensile;
        auto const& C_compressive = _ip_data[ip].C_compressive;

        auto const& strain_energy_tensile = _ip_data[ip].strain_energy_tensile;
        auto const& sigma_tensile = _ip_data[ip].sigma_tensile;

        auto& history_variable = _ip_data[ip].history_variable;
        auto& history_variable_prev = _ip_data[ip].history_variable_prev;

        auto const& sigma_real = _ip_data[ip].sigma_real;

        auto const k = _process_data.residual_stiffness(t, x_position)[0];
        auto const gc = _process_data.crack_resistance(t, x_position)[0];
        auto const ls = _process_data.crack_length_scale(t, x_position)[0];
        auto const M = _process_data.kinetic_coefficient(t, x_position)[0];
        auto rho_sr = _process_data.solid_density(t, x_position)[0];
        auto const alpha = _process_data.linear_thermal_expansion_coefficient(
            t, x_position)[0];
        double const c = _process_data.specific_heat_capacity(t, x_position)[0];
        auto const lambda =
            _process_data.thermal_conductivity(t, x_position)[0];
        auto const lambda_res =
            _process_data.residual_thermal_conductivity(t, x_position)[0];
        double const T0 = _process_data.reference_temperature;
        auto const& b = _process_data.specific_body_force;

        double const T_ip = N.dot(T);
        double const delta_T = T_ip - T0;
        // calculate real density
        double const rho_s = rho_sr * (1 - 3 * alpha * delta_T);

        // Kdd_component defines one term which both used in Kdd and local_rhs for phase
        // field
        typename ShapeMatricesType::NodalMatrixType const Kdd_component =
            dNdx.transpose() * 2 * gc * ls * dNdx;

        //
        // displacement equation, displacement part
        //

        double const d_ip = N.dot(d);

        double const degradation = d_ip * d_ip * (1 - k) + k;

        eps.noalias() = B * u;
        _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u, alpha,
                                                delta_T, degradation);

        local_Jac
            .template block<displacement_size, displacement_size>(
                displacement_index, displacement_index)
            .noalias() +=
            B.transpose() * (degradation * C_tensile + C_compressive) * B * w;

        typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                        displacement_size>
            N_u = ShapeMatricesType::template MatrixType<
                DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                          displacement_size);

        for (int i = 0; i < DisplacementDim; ++i)
            N_u.template block<1, displacement_size / DisplacementDim>(
                   i, i * displacement_size / DisplacementDim)
                .noalias() = N;

        using Invariants = MathLib::KelvinVector::Invariants<
                MathLib::KelvinVector::KelvinVectorDimensions<
                    DisplacementDim>::value>;

        local_rhs.template block<displacement_size, 1>(displacement_index, 0)
            .noalias() -=
            (B.transpose() * sigma_real - N_u.transpose() * rho_s * b) * w;

        //
        // displacement equation, temperature part
        //
        KuT.noalias() += B.transpose() *
                         (degradation * C_tensile + C_compressive) * alpha *
                         Invariants::identity2 * N * w;

        //
        // displacement equation, phasefield part
        //
        Kud.noalias() += B.transpose() * 2 * d_ip * sigma_tensile * N * w;

        if (history_variable_prev < strain_energy_tensile)
        {
            // INFO("History variable %g:", history_variable);
            // INFO("History variable previous %g:", history_variable_prev);
            history_variable = strain_energy_tensile;
            Kdu.noalias() +=
                N.transpose() * 2 * d_ip * sigma_tensile.transpose() * B * w;
            KdT.noalias() += N.transpose() * 2 * d_ip *
                             sigma_tensile.transpose() * alpha *
                             Invariants::identity2 * N * w;
        }
        else
        {
            history_variable = history_variable_prev;
        }

        //
        // phasefield equation, phasefield part
        //

        Kdd.noalias() += (Kdd_component + N.transpose() * 2 * history_variable * N +
                          N.transpose() * 0.5 * gc / ls * N) *
                         w;

        Ddd.noalias() += N.transpose() / M * N * w;

        //
        // temperature equation, temperature part;
        // temperature equation, phasefield part;
        // phasefield equation, temperature part
        //

        using Invariants = MathLib::KelvinVector::Invariants<
                MathLib::KelvinVector::KelvinVectorDimensions<
                    DisplacementDim>::value>;
        double const eps_m_trace = Invariants::trace(eps_m);
        if (eps_m_trace >= 0)
        {
            KTT.noalias() +=
                dNdx.transpose() *
                (d_ip * d_ip * lambda + (1 - d_ip) * (1 - d_ip) * lambda_res) *
                dNdx * w;
            KTd.noalias() +=
                dNdx.transpose() * 2 * d_ip * lambda * dNdx * T * N * w;
        }
        else
        {
            KTT.noalias() += dNdx.transpose() * lambda * dNdx * w;
        }

        DTT.noalias() += N.transpose() * rho_s * c * N * w;

        double const d_dot_ip = N.dot(d_dot);

        local_rhs.template block<phasefield_size, 1>(phasefield_index, 0)
            .noalias() -= (N.transpose() * d_dot_ip / M + Kdd_component * d +
                           N.transpose() * d_ip * 2 * history_variable -
                           N.transpose() * 0.5 * gc / ls * (1 - d_ip)) *
                          w;
    }
    // temperature equation, temperature part
    local_Jac
        .template block<temperature_size, temperature_size>(temperature_index,
                                                            temperature_index)
        .noalias() += KTT + DTT / dt;
    // temperature equation, phasefield part
    local_Jac
        .template block<temperature_size, phasefield_size>(temperature_index,
                                                           phasefield_index)
        .noalias() += KTd;
    // displacement equation, temperature part
    local_Jac
        .template block<displacement_size, temperature_size>(displacement_index,
                                                             temperature_index)
        .noalias() -= KuT;
    // displacement equation, phasefield part
    local_Jac
        .template block<displacement_size, phasefield_size>(displacement_index,
                                                            phasefield_index)
        .noalias() += Kud;
    // phasefield equation, phasefield part.
    local_Jac
        .template block<phasefield_size, phasefield_size>(phasefield_index,
                                                          phasefield_index)
        .noalias() += Kdd + Ddd / dt;
    // phasefield equation, displacement part.
    local_Jac
        .template block<phasefield_size, displacement_size>(phasefield_index,
                                                            displacement_index)
        .noalias() += Kdu;
    // phasefield equation, temperature part
    local_Jac
        .template block<phasefield_size, temperature_size>(phasefield_index,
                                                           temperature_index)
        .noalias() -= KdT;

    local_rhs.template block<temperature_size, 1>(temperature_index, 0)
        .noalias() -= KTT * T + DTT * T_dot;

    local_rhs.template block<temperature_size, 1>(temperature_index, 0)
        .noalias() -= KTd * d;

    local_rhs.template block<phasefield_size, 1>(phasefield_index, 0)
        .noalias() += KdT * T;

    local_rhs.template block<displacement_size, 1>(displacement_index, 0)
        .noalias() += KuT * T;
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void ThermoMechanicalPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                DisplacementDim>::
    assembleWithJacobianForStaggeredScheme(
        double const t, std::vector<double> const& local_xdot,
        const double dxdot_dx, const double dx_dx,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions)
{
    if (local_coupled_solutions.process_id == _phase_field_process_id)
    {
        assembleWithJacobianForPhaseFieldEquations(
            t, local_xdot, dxdot_dx, dx_dx, local_M_data, local_K_data,
            local_b_data, local_Jac_data, local_coupled_solutions);
        return;
    }

    if (local_coupled_solutions.process_id == _heat_conduction_process_id)
    {
        assembleWithJacobianForHeatConductionEquations(
            t, local_xdot, dxdot_dx, dx_dx, local_M_data, local_K_data,
            local_b_data, local_Jac_data, local_coupled_solutions);
        return;
    }

    // For the equations with deformation
    assembleWithJacobianForDeformationEquations(
        t, local_xdot, dxdot_dx, dx_dx, local_M_data, local_K_data,
        local_b_data, local_Jac_data, local_coupled_solutions);
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void ThermoMechanicalPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                DisplacementDim>::
    assembleWithJacobianForDeformationEquations(
        double const t, std::vector<double> const& /*local_xdot*/,
        const double /*dxdot_dx*/, const double /*dx_dx*/,
        std::vector<double>& /*local_M_data*/, std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions)
{
    auto const& local_d =
        local_coupled_solutions.local_coupled_xs[_phase_field_process_id];
    auto const& local_u =
        local_coupled_solutions.local_coupled_xs[_mechanics_related_process_id];
    auto const& local_T =
        local_coupled_solutions.local_coupled_xs[_heat_conduction_process_id];
    assert(local_T.size() == temperature_size);
    assert(local_d.size() == phasefield_size);
    assert(local_u.size() == displacement_size);

    auto d = Eigen::Map<
        typename ShapeMatricesType::template VectorType<phasefield_size> const>(
        local_d.data(), phasefield_size);

    auto u = Eigen::Map<typename ShapeMatricesType::template VectorType<
        displacement_size> const>(local_u.data(), displacement_size);

    auto T = Eigen::Map<typename ShapeMatricesType::template VectorType<
        temperature_size> const>(local_T.data(), temperature_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesType::template MatrixType<
            displacement_size, displacement_size>>(
        local_Jac_data, displacement_size, displacement_size);

    auto local_rhs =
        MathLib::createZeroedVector<typename ShapeMatricesType::
                                        template VectorType<displacement_size>>(
            local_b_data, displacement_size);

    double const& dt = _process_data.dt;

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        auto const x_coord =
            interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(_element,
                                                                     N);
        auto const& B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunction::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx, N, x_coord, _is_axially_symmetric);

        auto& eps = _ip_data[ip].eps;

        auto const& C_tensile = _ip_data[ip].C_tensile;
        auto const& C_compressive = _ip_data[ip].C_compressive;

        auto const& sigma_real = _ip_data[ip].sigma_real;

        auto const k = _process_data.residual_stiffness(t, x_position)[0];
        auto rho_sr = _process_data.solid_density(t, x_position)[0];
        auto const alpha = _process_data.linear_thermal_expansion_coefficient(
            t, x_position)[0];
        double const T0 = _process_data.reference_temperature;
        auto const& b = _process_data.specific_body_force;

        double const T_ip = N.dot(T);
        double const delta_T = T_ip - T0;
        // calculate real density
        double const rho_s = rho_sr * (1 - 3 * alpha * delta_T);

        double const d_ip = N.dot(d);
        double const degradation = d_ip * d_ip * (1 - k) + k;

        auto const C_eff = degradation * C_tensile + C_compressive;
        eps.noalias() = B * u;
        _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u, alpha,
                                                delta_T, degradation);

        typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                        displacement_size>
            N_u = ShapeMatricesType::template MatrixType<
                DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                          displacement_size);

        for (int i = 0; i < DisplacementDim; ++i)
            N_u.template block<1, displacement_size / DisplacementDim>(
                   i, i * displacement_size / DisplacementDim)
                .noalias() = N;

        using Invariants = MathLib::KelvinVector::Invariants<
                MathLib::KelvinVector::KelvinVectorDimensions<
                    DisplacementDim>::value>;

        local_Jac.noalias() +=
            B.transpose() * C_eff * B * w;

        local_rhs.noalias() -=
            (B.transpose() * sigma_real - N_u.transpose() * rho_s * b) * w;

    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void ThermoMechanicalPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                DisplacementDim>::
    assembleWithJacobianForHeatConductionEquations(
        double const t, std::vector<double> const& local_xdot,
        const double /*dxdot_dx*/, const double /*dx_dx*/,
        std::vector<double>& /*local_M_data*/, std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions)
{
    auto const& local_d =
        local_coupled_solutions.local_coupled_xs[_phase_field_process_id];
    auto const& local_T =
        local_coupled_solutions.local_coupled_xs[_heat_conduction_process_id];
    assert(local_T.size() == temperature_size);
    assert(local_d.size() == phasefield_size);

    auto d = Eigen::Map<
        typename ShapeMatricesType::template VectorType<phasefield_size> const>(
        local_d.data(), phasefield_size);

    auto T = Eigen::Map<typename ShapeMatricesType::template VectorType<
        temperature_size> const>(local_T.data(), temperature_size);

    auto T_dot = Eigen::Map<typename ShapeMatricesType::template VectorType<
        temperature_size> const>(local_xdot.data(), temperature_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesType::template MatrixType<
            temperature_size, temperature_size>>(
        local_Jac_data, temperature_size, temperature_size);

    auto local_rhs =
        MathLib::createZeroedVector<typename ShapeMatricesType::
                                        template VectorType<temperature_size>>(
            local_b_data, temperature_size);

    double const& dt = _process_data.dt;

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        auto& eps_m = _ip_data[ip].eps_m;
        auto& heatflux = _ip_data[ip].heatflux;

        auto rho_sr = _process_data.solid_density(t, x_position)[0];
        auto const alpha = _process_data.linear_thermal_expansion_coefficient(
            t, x_position)[0];
        double const c = _process_data.specific_heat_capacity(t, x_position)[0];
        auto const lambda =
            _process_data.thermal_conductivity(t, x_position)[0];
        auto const lambda_res =
            _process_data.residual_thermal_conductivity(t, x_position)[0];
        double const T0 = _process_data.reference_temperature;

        double const d_ip = N.dot(d);
        double const T_ip = N.dot(T);
        double const T_dot_ip = N.dot(T_dot);
        double const delta_T = T_ip - T0;
        // calculate real density
        double const rho_s = rho_sr * (1 - 3 * alpha * delta_T);
        // calculate effective thermal conductivity
        auto const lambda_eff = d_ip * d_ip * lambda +
                                (1 - d_ip) * (1 - d_ip) * lambda_res;

        using Invariants = MathLib::KelvinVector::Invariants<
                MathLib::KelvinVector::KelvinVectorDimensions<
                    DisplacementDim>::value>;

        double const eps_m_trace = Invariants::trace(eps_m);
        if (eps_m_trace >= 0)
        {
            local_Jac.noalias() +=
                (dNdx.transpose() * lambda_eff * dNdx +
                 N.transpose() * rho_s * c * N / dt) * w;

            local_rhs.noalias() -= (N.transpose() * rho_s * c * T_dot_ip +
                    dNdx.transpose() * lambda_eff * dNdx * T) * w;

            heatflux.noalias() = - (lambda_eff * dNdx * T) * w;
        }
        else
        {
            local_Jac.noalias() += (dNdx.transpose() * lambda * dNdx +
                                    N.transpose() * rho_s * c * N / dt) * w;

            local_rhs.noalias() -= (N.transpose() * rho_s * c * T_dot_ip +
                    dNdx.transpose() * lambda * dNdx * T) * w;

            heatflux.noalias() = - (lambda * dNdx * T) * w;
        }
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void ThermoMechanicalPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                DisplacementDim>::
    assembleWithJacobianForPhaseFieldEquations(
        double const t, std::vector<double> const& /*local_xdot*/,
        const double /*dxdot_dx*/, const double /*dx_dx*/,
        std::vector<double>& /*local_M_data*/, std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions)
{
    auto const& local_d =
        local_coupled_solutions.local_coupled_xs[_phase_field_process_id];
    assert(local_d.size() == phasefield_size);

    auto d = Eigen::Map<
        typename ShapeMatricesType::template VectorType<phasefield_size> const>(
        local_d.data(), phasefield_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesType::template MatrixType<
            phasefield_size, phasefield_size>>(
        local_Jac_data, phasefield_size, phasefield_size);

    auto local_rhs =
        MathLib::createZeroedVector<typename ShapeMatricesType::
                                        template VectorType<phasefield_size>>(
            local_b_data, phasefield_size);

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        auto const& strain_energy_tensile = _ip_data[ip].strain_energy_tensile;

        auto const gc = _process_data.crack_resistance(t, x_position)[0];
        auto const ls = _process_data.crack_length_scale(t, x_position)[0];

        double const d_ip = N.dot(d);

        using Invariants = MathLib::KelvinVector::Invariants<
                MathLib::KelvinVector::KelvinVectorDimensions<
                    DisplacementDim>::value>;

        local_Jac.noalias() +=
                (dNdx.transpose() * gc * ls * dNdx +
                 N.transpose() * 2 * strain_energy_tensile * N +
                 N.transpose() * gc / ls * N) * w;

        local_rhs.noalias() -= (dNdx.transpose() * gc * ls * dNdx * d +
                        N.transpose() * d_ip * 2 * strain_energy_tensile -
                        N.transpose() * gc / ls * (1 - d_ip)) * w;

    }
}

}  // namespace ThermoMechanicalPhaseField
}  // namespace ProcessLib
