/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file ThermoMechanicsFEM-impl.h
 *
 * Created on July 2, 2019, 2:12 PM
 *
 */

#pragma once

#include "ThermoMechanicsFEM.h"

namespace ProcessLib
{
namespace ThermoMechanics
{
template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void ThermoMechanicsLocalAssembler<ShapeFunction, IntegrationMethod,
                                   DisplacementDim>::
    assembleWithJacobianForStaggeredScheme(
        const double t,
        const std::vector<double>& local_xdot,
        const double dxdot_dx,
        const double dx_dx,
        std::vector<double>& local_M_data,
        std::vector<double>& local_K_data,
        std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data,
        const LocalCoupledSolutions& local_coupled_solutions)
{
    // For the equations with pressure
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
void ThermoMechanicsLocalAssembler<ShapeFunction, IntegrationMethod,
                                   DisplacementDim>::
    assembleWithJacobianForDeformationEquations(
        const double t, const std::vector<double>& /*local_xdot*/,
        const double /*dxdot_dx*/, const double /*dx_dx*/,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        const LocalCoupledSolutions& local_coupled_solutions)
{
    auto const& local_T_vector =
        local_coupled_solutions.local_coupled_xs[_heat_conduction_process_id];
    assert(local_T_vector.size() == temperature_size);
    auto const local_T =
        Eigen::Map<typename ShapeMatricesType::template VectorType<
            temperature_size> const>(local_T_vector.data(), temperature_size);

    auto const& local_T0_vector = local_coupled_solutions.local_coupled_xs0[0];
    assert(local_T0_vector.size() == temperature_size);
    auto const local_T0 =
        Eigen::Map<typename ShapeMatricesType::template VectorType<
            temperature_size> const>(local_T0_vector.data(), temperature_size);

    auto const& local_u_vector =
        local_coupled_solutions.local_coupled_xs[_mechanics_process_id];
    assert(local_u_vector.size() == displacement_size);
    auto const u = Eigen::Map<typename ShapeMatricesType::template VectorType<
        displacement_size> const>(local_u_vector.data(), displacement_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesType::template MatrixType<displacement_size,
                                                        displacement_size>>(
        local_Jac_data, displacement_size, displacement_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesType::template VectorType<displacement_size>>(
        local_b_data, displacement_size);

    double const& dt = _process_data.dt;

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
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

        auto& sigma = _ip_data[ip].sigma;
        auto const& sigma_prev = _ip_data[ip].sigma_prev;

        auto& eps = _ip_data[ip].eps;
        auto const& eps_prev = _ip_data[ip].eps_prev;

        auto& eps_m = _ip_data[ip].eps_m;
        auto const& eps_m_prev = _ip_data[ip].eps_m_prev;

        auto& state = _ip_data[ip].material_state_variables;

        const double T = N.dot(local_T);  // T at integration point
        double const dT = T - N.dot(local_T0);
        // calculate thermally induced strain
        // assume isotropic thermal expansion
        auto const alpha = _process_data.linear_thermal_expansion_coefficient(
            t, x_position)[0];
        double const linear_thermal_strain_increment = alpha * dT;

        //
        // displacement equation, displacement part
        //
        eps.noalias() = B * u;

        using Invariants = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value>;

        // assume isotropic thermal expansion
        eps_m.noalias() =
            eps_m_prev + eps - eps_prev -
            linear_thermal_strain_increment * Invariants::identity2;

        auto&& solution = _ip_data[ip].solid_material.integrateStress(
            t, x_position, dt, eps_m_prev, eps_m, sigma_prev, *state, T);

        if (!solution)
        {
            OGS_FATAL("Computation of local constitutive relation failed.");
        }

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
        std::tie(sigma, state, C) = std::move(*solution);

        local_Jac.noalias() += B.transpose() * C * B * w;

        typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                        displacement_size>
            N_u = ShapeMatricesType::template MatrixType<
                DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                          displacement_size);

        for (int i = 0; i < DisplacementDim; ++i)
        {
            N_u.template block<1, displacement_size / DisplacementDim>(
                   i, i * displacement_size / DisplacementDim)
                .noalias() = N;
        }

        // calculate real density
        // rho_s_{n+1} * (V_{n} + dV) = rho_s_n * V_n
        // dV = 3 * alpha * dT * V_0
        // rho_s_{n+1} = rho_s_n / (1 + 3 * alpha * dT )
        // see reference solid density description for details.
        auto& rho_s = _ip_data[ip].solid_density;
        rho_s = _ip_data[ip].solid_density_prev /
                (1 + 3 * linear_thermal_strain_increment);

        auto const& b = _process_data.specific_body_force;
        local_rhs.noalias() -=
            (B.transpose() * sigma - N_u.transpose() * rho_s * b) * w;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void ThermoMechanicsLocalAssembler<ShapeFunction, IntegrationMethod,
                                   DisplacementDim>::
    assembleWithJacobianForHeatConductionEquations(
        const double t, const std::vector<double>& /*local_xdot*/,
        const double /*dxdot_dx*/, const double /*dx_dx*/,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        const LocalCoupledSolutions& local_coupled_solutions)
{
    auto const& local_T_vector =
        local_coupled_solutions.local_coupled_xs[_heat_conduction_process_id];
    assert(local_T_vector.size() == temperature_size);
    auto const local_T =
        Eigen::Map<typename ShapeMatricesType::template VectorType<
            temperature_size> const>(local_T_vector.data(), temperature_size);

    auto const& local_T0_vector = local_coupled_solutions.local_coupled_xs0[0];
    assert(local_T0_vector.size() == temperature_size);
    auto const local_dT =
        local_T -
        Eigen::Map<typename ShapeMatricesType::template VectorType<
            temperature_size> const>(local_T0_vector.data(), temperature_size);

    auto local_Jac =
        MathLib::createZeroedMatrix <
        typename ShapeMatricesType::template MatrixType<temperature_size,
                                                        temperature_size>>(
            local_Jac_data, temperature_size, temperature_size);

    auto local_rhs =
        MathLib::createZeroedVector <
        typename ShapeMatricesType::template VectorType<temperature_size>>(
            local_b_data, temperature_size);

    typename ShapeMatricesType::NodalMatrixType mass;
    mass.setZero(temperature_size, temperature_size);

    typename ShapeMatricesType::NodalMatrixType laplace;
    laplace.setZero(temperature_size, temperature_size);

    double const& dt = _process_data.dt;

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        // calculate real density
        // rho_s_{n+1} * (V_{n} + dV) = rho_s_n * V_n
        // dV = 3 * alpha * dT * V_0
        // rho_s_{n+1} = rho_s_n / (1 + 3 * alpha * dT )
        // see reference solid density description for details.
        auto& rho_s = _ip_data[ip].solid_density;
        // calculate thermally induced strain
        // assume isotropic thermal expansion
        auto const alpha = _process_data.linear_thermal_expansion_coefficient(
            t, x_position)[0];

        double const dT = N.dot(local_dT);
        double const linear_thermal_strain_increment = alpha * dT;
        rho_s = _ip_data[ip].solid_density_prev /
                (1 + 3 * linear_thermal_strain_increment);
        auto const c_p = _process_data.specific_heat_capacity(t, x_position)[0];
        mass.noalias() += N.transpose() * rho_s * c_p * N * w;

        auto const lambda =
            _process_data.thermal_conductivity(t, x_position)[0];
        laplace.noalias() += dNdx.transpose() * lambda * dNdx * w;
    }
    local_Jac.noalias() += laplace + mass / dt;

    local_rhs.noalias() -= laplace * local_T + mass * local_dT / dt;
}

}  // namespace ThermoMechanics
}  // namespace ProcessLib
