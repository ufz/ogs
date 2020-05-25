/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
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
    assembleWithJacobianForStaggeredScheme(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        std::vector<double> const& local_xdot, const double dxdot_dx,
        const double dx_dx, int const process_id,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& /*local_coupled_solutions*/)
{
    if (process_id == phase_field_process_id_)
    {
        assembleWithJacobianForPhaseFieldEquations(
            t, local_x, local_xdot, dxdot_dx, dx_dx, local_M_data, local_K_data,
            local_b_data, local_Jac_data);
        return;
    }

    if (process_id == heat_conduction_process_id_)
    {
        assembleWithJacobianForHeatConductionEquations(
            t, dt, local_x, local_xdot, dxdot_dx, dx_dx, local_M_data,
            local_K_data, local_b_data, local_Jac_data);
        return;
    }

    // For the equations with deformation
    assembleWithJacobianForDeformationEquations(
        t, dt, local_x, local_xdot, dxdot_dx, dx_dx, local_M_data, local_K_data,
        local_b_data, local_Jac_data);
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void ThermoMechanicalPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                              DisplacementDim>::
    assembleWithJacobianForDeformationEquations(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        std::vector<double> const& /*local_xdot*/, const double /*dxdot_dx*/,
        const double /*dx_dx*/, std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    assert(local_x.size() ==
           phasefield_size + displacement_size + temperature_size);

    auto const d = local_x.template segment<phasefield_size>(phasefield_index);
    auto const u =
        local_x.template segment<displacement_size>(displacement_index);
    auto const T =
        local_x.template segment<temperature_size>(temperature_index);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesType::template MatrixType<displacement_size,
                                                        displacement_size>>(
        local_Jac_data, displacement_size, displacement_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesType::template VectorType<displacement_size>>(
        local_b_data, displacement_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(element_.getID());

    int const n_integration_points = integration_method_.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = ip_data_[ip].integration_weight;
        auto const& N = ip_data_[ip].N;
        auto const& dNdx = ip_data_[ip].dNdx;

        auto const x_coord =
            interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(element_,
                                                                     N);
        auto const& B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunction::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx, N, x_coord, is_axially_symmetric_);

        auto& eps = ip_data_[ip].eps;

        auto const& C_tensile = ip_data_[ip].C_tensile;
        auto const& C_compressive = ip_data_[ip].C_compressive;

        auto const& sigma = ip_data_[ip].sigma;

        auto const k = process_data_.residual_stiffness(t, x_position)[0];
        auto rho_sr = process_data_.solid_density(t, x_position)[0];
        auto const alpha = process_data_.linear_thermal_expansion_coefficient(
            t, x_position)[0];
        double const T0 = process_data_.reference_temperature;
        auto const& b = process_data_.specific_body_force;

        double const T_ip = N.dot(T);
        double const delta_T = T_ip - T0;
        // calculate real density
        double const rho_s = rho_sr / (1 + 3 * alpha * delta_T);

        double const d_ip = N.dot(d);
        double const degradation = d_ip * d_ip * (1 - k) + k;

        auto const C_eff = degradation * C_tensile + C_compressive;
        eps.noalias() = B * u;
        ip_data_[ip].updateConstitutiveRelation(t, x_position, dt, u, alpha,
                                                delta_T, degradation);

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

        local_Jac.noalias() += B.transpose() * C_eff * B * w;

        local_rhs.noalias() -=
            (B.transpose() * sigma - N_u.transpose() * rho_s * b) * w;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void ThermoMechanicalPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                              DisplacementDim>::
    assembleWithJacobianForHeatConductionEquations(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        std::vector<double> const& local_xdot, const double /*dxdot_dx*/,
        const double /*dx_dx*/, std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    assert(local_x.size() ==
           phasefield_size + displacement_size + temperature_size);

    auto const d = local_x.template segment<phasefield_size>(phasefield_index);
    auto const T =
        local_x.template segment<temperature_size>(temperature_index);

    auto T_dot = Eigen::Map<typename ShapeMatricesType::template VectorType<
        temperature_size> const>(local_xdot.data(), temperature_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesType::template MatrixType<temperature_size,
                                                        temperature_size>>(
        local_Jac_data, temperature_size, temperature_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesType::template VectorType<temperature_size>>(
        local_b_data, temperature_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(element_.getID());

    int const n_integration_points = integration_method_.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = ip_data_[ip].integration_weight;
        auto const& N = ip_data_[ip].N;
        auto const& dNdx = ip_data_[ip].dNdx;

        auto& eps_m = ip_data_[ip].eps_m;
        auto& heatflux = ip_data_[ip].heatflux;

        auto rho_sr = process_data_.solid_density(t, x_position)[0];
        auto const alpha = process_data_.linear_thermal_expansion_coefficient(
            t, x_position)[0];
        double const c = process_data_.specific_heat_capacity(t, x_position)[0];
        auto const lambda =
            process_data_.thermal_conductivity(t, x_position)[0];
        auto const lambda_res =
            process_data_.residual_thermal_conductivity(t, x_position)[0];
        double const T0 = process_data_.reference_temperature;

        double const d_ip = N.dot(d);
        double const T_ip = N.dot(T);
        double const T_dot_ip = N.dot(T_dot);
        double const delta_T = T_ip - T0;
        // calculate real density
        double const rho_s = rho_sr / (1 + 3 * alpha * delta_T);
        // calculate effective thermal conductivity
        auto const lambda_eff =
            d_ip * d_ip * lambda + (1 - d_ip) * (1 - d_ip) * lambda_res;

        using Invariants = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value>;

        double const eps_m_trace = Invariants::trace(eps_m);
        if (eps_m_trace >= 0)
        {
            local_Jac.noalias() += (dNdx.transpose() * lambda_eff * dNdx +
                                    N.transpose() * rho_s * c * N / dt) *
                                   w;

            local_rhs.noalias() -= (N.transpose() * rho_s * c * T_dot_ip +
                                    dNdx.transpose() * lambda_eff * dNdx * T) *
                                   w;

            heatflux.noalias() = -(lambda_eff * dNdx * T) * w;
        }
        else
        {
            local_Jac.noalias() += (dNdx.transpose() * lambda * dNdx +
                                    N.transpose() * rho_s * c * N / dt) *
                                   w;

            local_rhs.noalias() -= (N.transpose() * rho_s * c * T_dot_ip +
                                    dNdx.transpose() * lambda * dNdx * T) *
                                   w;

            heatflux.noalias() = -(lambda * dNdx * T) * w;
        }
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void ThermoMechanicalPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                              DisplacementDim>::
    assembleWithJacobianForPhaseFieldEquations(
        double const t, Eigen::VectorXd const& local_x,
        std::vector<double> const& /*local_xdot*/, const double /*dxdot_dx*/,
        const double /*dx_dx*/, std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    assert(local_x.size() ==
           phasefield_size + displacement_size + temperature_size);

    auto const d = local_x.template segment<phasefield_size>(phasefield_index);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesType::template MatrixType<phasefield_size,
                                                        phasefield_size>>(
        local_Jac_data, phasefield_size, phasefield_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesType::template VectorType<phasefield_size>>(
        local_b_data, phasefield_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(element_.getID());

    int const n_integration_points = integration_method_.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = ip_data_[ip].integration_weight;
        auto const& N = ip_data_[ip].N;
        auto const& dNdx = ip_data_[ip].dNdx;

        auto const& strain_energy_tensile = ip_data_[ip].strain_energy_tensile;

        auto const gc = process_data_.crack_resistance(t, x_position)[0];
        auto const ls = process_data_.crack_length_scale(t, x_position)[0];

        double const d_ip = N.dot(d);

        local_Jac.noalias() += (dNdx.transpose() * gc * ls * dNdx +
                                N.transpose() * 2 * strain_energy_tensile * N +
                                N.transpose() * gc / ls * N) *
                               w;

        local_rhs.noalias() -=
            (dNdx.transpose() * gc * ls * dNdx * d +
             N.transpose() * d_ip * 2 * strain_energy_tensile -
             N.transpose() * gc / ls * (1 - d_ip)) *
            w;
    }
}

}  // namespace ThermoMechanicalPhaseField
}  // namespace ProcessLib
