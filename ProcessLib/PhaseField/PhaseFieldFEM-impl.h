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

#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
namespace PhaseField
{
template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void PhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                              DisplacementDim>::
    assembleWithJacobianForStaggeredScheme(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        std::vector<double> const& local_xdot, const double dxdot_dx,
        const double dx_dx, int const process_id,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& /*local_coupled_solutions*/)
{
    // For the equations with phase field.
    if (process_id == 1)
    {
        assembleWithJacobianPhaseFieldEquations(
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
void PhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                              DisplacementDim>::
    assembleWithJacobianForDeformationEquations(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        std::vector<double> const& /*local_xdot*/, const double /*dxdot_dx*/,
        const double /*dx_dx*/, std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    using DeformationMatrix =
        typename ShapeMatricesType::template MatrixType<displacement_size,
                                                        displacement_size>;

    auto const d = local_x.template segment<phasefield_size>(phasefield_index);
    auto const u =
        local_x.template segment<displacement_size>(displacement_index);

    auto local_Jac = MathLib::createZeroedMatrix<DeformationMatrix>(
        local_Jac_data, displacement_size, displacement_size);

    auto local_rhs = MathLib::createZeroedVector<DeformationVector>(
        local_b_data, displacement_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(element_.getID());

    auto local_pressure = 0.0;
    if (process_data_.crack_pressure)
    {
        local_pressure = process_data_.unity_pressure;
    }

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
        eps.noalias() = B * u;
        double const k = process_data_.residual_stiffness(t, x_position)[0];
        double const d_ip = N.dot(d);
        double const degradation = d_ip * d_ip * (1 - k) + k;
        ip_data_[ip].updateConstitutiveRelation(t, x_position, dt, u,
                                                degradation);

        auto& sigma = ip_data_[ip].sigma;
        auto const& C_tensile = ip_data_[ip].C_tensile;
        auto const& C_compressive = ip_data_[ip].C_compressive;

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

        auto const rho_sr = process_data_.solid_density(t, x_position)[0];
        auto const& b = process_data_.specific_body_force;

        local_rhs.noalias() -=
            (B.transpose() * sigma - N_u.transpose() * rho_sr * b -
             local_pressure * N_u.transpose() * dNdx * d) *
            w;

        local_Jac.noalias() +=
            B.transpose() * (degradation * C_tensile + C_compressive) * B * w;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void PhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                              DisplacementDim>::
    assembleWithJacobianPhaseFieldEquations(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        std::vector<double> const& /*local_xdot*/, const double /*dxdot_dx*/,
        const double /*dx_dx*/, std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    auto const d = local_x.template segment<phasefield_size>(phasefield_index);
    auto const u =
        local_x.template segment<displacement_size>(displacement_index);

    auto local_Jac = MathLib::createZeroedMatrix<PhaseFieldMatrix>(
        local_Jac_data, phasefield_size, phasefield_size);
    auto local_rhs = MathLib::createZeroedVector<PhaseFieldVector>(
        local_b_data, phasefield_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(element_.getID());

    auto local_pressure = 0.0;
    if (process_data_.crack_pressure)
    {
        local_pressure = process_data_.unity_pressure;
    }
    else if (process_data_.propagating_crack)
    {
        local_pressure = process_data_.pressure;
    }

    int const n_integration_points = integration_method_.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = ip_data_[ip].integration_weight;
        auto const& N = ip_data_[ip].N;
        auto const& dNdx = ip_data_[ip].dNdx;

        double const gc = process_data_.crack_resistance(t, x_position)[0];
        double const ls = process_data_.crack_length_scale(t, x_position)[0];

        // for propagating crack, u is rescaled.
        if (process_data_.propagating_crack)
        {
            double const k = process_data_.residual_stiffness(t, x_position)[0];
            double const d_ip = N.dot(d);
            double const degradation = d_ip * d_ip * (1 - k) + k;
            auto const x_coord =
                interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                    element_, N);
            auto const& B = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunction::NPOINTS,
                typename BMatricesType::BMatrixType>(dNdx, N, x_coord,
                                                     is_axially_symmetric_);

            auto& eps = ip_data_[ip].eps;
            eps.noalias() = B * u;
            ip_data_[ip].updateConstitutiveRelation(t, x_position, dt, u,
                                                    degradation);
        }

        auto const& strain_energy_tensile = ip_data_[ip].strain_energy_tensile;

        auto& ip_data = ip_data_[ip];
        ip_data.strain_energy_tensile = strain_energy_tensile;

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

        local_Jac.noalias() +=
            (2 * N.transpose() * N * strain_energy_tensile +
             gc * (N.transpose() * N / ls + dNdx.transpose() * dNdx * ls)) *
            w;

        local_rhs.noalias() -=
            (N.transpose() * N * d * 2 * strain_energy_tensile +
             gc * ((N.transpose() * N / ls + dNdx.transpose() * dNdx * ls) * d -
                   N.transpose() / ls) -
             local_pressure * dNdx.transpose() * N_u * u) *
            w;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void PhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                              DisplacementDim>::
    computeCrackIntegral(std::size_t mesh_item_id,
                         std::vector<std::reference_wrapper<
                             NumLib::LocalToGlobalIndexMap>> const& dof_tables,
                         GlobalVector const& /*x*/, double const /*t*/,
                         double& crack_volume,
                         CoupledSolutionsForStaggeredScheme const* const cpl_xs)
{
    assert(cpl_xs != nullptr);

    std::vector<std::vector<GlobalIndexType>> indices_of_processes;
    indices_of_processes.reserve(dof_tables.size());
    std::transform(dof_tables.begin(), dof_tables.end(),
                   std::back_inserter(indices_of_processes),
                   [&](NumLib::LocalToGlobalIndexMap const& dof_table) {
                       return NumLib::getIndices(mesh_item_id, dof_table);
                   });

    auto local_coupled_xs =
        getCoupledLocalSolutions(cpl_xs->coupled_xs, indices_of_processes);
    assert(local_coupled_xs.size() == displacement_size + phasefield_size);

    auto const d = Eigen::Map<PhaseFieldVector const>(
        &local_coupled_xs[phasefield_index], phasefield_size);
    auto const u = Eigen::Map<DeformationVector const>(
        &local_coupled_xs[displacement_index], displacement_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(element_.getID());

    int const n_integration_points = integration_method_.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = ip_data_[ip].integration_weight;
        auto const& N = ip_data_[ip].N;
        auto const& dNdx = ip_data_[ip].dNdx;

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

        crack_volume += (N_u * u).dot(dNdx * d) * w;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void PhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                              DisplacementDim>::
    computeEnergy(std::size_t mesh_item_id,
                  std::vector<std::reference_wrapper<
                      NumLib::LocalToGlobalIndexMap>> const& dof_tables,
                  GlobalVector const& /*x*/, double const t,
                  double& elastic_energy, double& surface_energy,
                  double& pressure_work,
                  CoupledSolutionsForStaggeredScheme const* const cpl_xs)
{
    assert(cpl_xs != nullptr);

    std::vector<std::vector<GlobalIndexType>> indices_of_processes;
    indices_of_processes.reserve(dof_tables.size());
    std::transform(dof_tables.begin(), dof_tables.end(),
                   std::back_inserter(indices_of_processes),
                   [&](NumLib::LocalToGlobalIndexMap const& dof_table) {
                       return NumLib::getIndices(mesh_item_id, dof_table);
                   });

    auto const local_coupled_xs =
        getCoupledLocalSolutions(cpl_xs->coupled_xs, indices_of_processes);
    assert(local_coupled_xs.size() == displacement_size + phasefield_size);

    auto const d = Eigen::Map<PhaseFieldVector const>(
        &local_coupled_xs[phasefield_index], phasefield_size);
    auto const u = Eigen::Map<DeformationVector const>(
        &local_coupled_xs[displacement_index], displacement_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(element_.getID());

    int const n_integration_points = integration_method_.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = ip_data_[ip].integration_weight;
        auto const& N = ip_data_[ip].N;
        auto const& dNdx = ip_data_[ip].dNdx;
        double const d_ip = N.dot(d);
        auto pressure_ip = process_data_.pressure;
        auto u_corrected = pressure_ip * u;
        double const gc = process_data_.crack_resistance(t, x_position)[0];
        double const ls = process_data_.crack_length_scale(t, x_position)[0];

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

        elastic_energy += ip_data_[ip].elastic_energy * w;

        surface_energy +=
            0.5 * gc *
            ((1 - d_ip) * (1 - d_ip) / ls + (dNdx * d).dot((dNdx * d)) * ls) *
            w;

        if (process_data_.crack_pressure)
        {
            pressure_work +=
                pressure_ip * (N_u * u_corrected).dot(dNdx * d) * w;
        }
    }
}
}  // namespace PhaseField
}  // namespace ProcessLib
