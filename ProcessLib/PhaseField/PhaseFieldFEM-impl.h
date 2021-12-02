/**
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
        Eigen::VectorXd const& /*local_xdot*/, int const process_id,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    // For the equations with phase field.
    if (process_id == phase_process_id)
    {
        assembleWithJacobianPhaseFieldEquations(t, dt, local_x, local_b_data,
                                                local_Jac_data);
        return;
    }

    // For the equations with deformation
    assembleWithJacobianForDeformationEquations(t, dt, local_x, local_b_data,
                                                local_Jac_data);
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void PhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                              DisplacementDim>::
    assembleWithJacobianForDeformationEquations(
        double const t, double const dt, Eigen::VectorXd const& local_x,
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
    x_position.setElementID(_element.getID());

    auto local_pressure = 0.0;
    if (_process_data.crack_pressure)
    {
        local_pressure = _process_data.unity_pressure;
    }

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                _element, N);

        auto const& B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunction::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx, N, x_coord, _is_axially_symmetric);

        auto& eps = _ip_data[ip].eps;
        eps.noalias() = B * u;
        double const k = _process_data.residual_stiffness(t, x_position)[0];
        double const d_ip = N.dot(d);
        double const degradation = d_ip * d_ip * (1 - k) + k;
        _ip_data[ip].updateConstitutiveRelation(
            t, x_position, dt, u, degradation,
            _process_data.energy_split_model);

        auto& sigma = _ip_data[ip].sigma;
        auto const& C_tensile = _ip_data[ip].C_tensile;
        auto const& C_compressive = _ip_data[ip].C_compressive;

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

        auto const rho_sr = _process_data.solid_density(t, x_position)[0];
        auto const& b = _process_data.specific_body_force;

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
    assembleWithJacobianPhaseFieldEquations(double const t, double const dt,
                                            Eigen::VectorXd const& local_x,
                                            std::vector<double>& local_b_data,
                                            std::vector<double>& local_Jac_data)
{
    auto const d = local_x.template segment<phasefield_size>(phasefield_index);
    auto const u =
        local_x.template segment<displacement_size>(displacement_index);

    auto local_Jac = MathLib::createZeroedMatrix<PhaseFieldMatrix>(
        local_Jac_data, phasefield_size, phasefield_size);
    auto local_rhs = MathLib::createZeroedVector<PhaseFieldVector>(
        local_b_data, phasefield_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    auto local_pressure = 0.0;
    if (_process_data.crack_pressure)
    {
        local_pressure = _process_data.unity_pressure;
    }
    else if (_process_data.hydro_crack)
    {
        local_pressure = _process_data.pressure;
    }

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        double const gc = _process_data.crack_resistance(t, x_position)[0];
        double const ls = _process_data.crack_length_scale(t, x_position)[0];

        double const k = _process_data.residual_stiffness(t, x_position)[0];
        double const d_ip = N.dot(d);
        double const degradation = d_ip * d_ip * (1 - k) + k;
        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                _element, N);
        auto const& B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunction::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx, N, x_coord, _is_axially_symmetric);

        auto& eps = _ip_data[ip].eps;
        eps.noalias() = B * u;
        _ip_data[ip].updateConstitutiveRelation(
            t, x_position, dt, u, degradation,
            _process_data.energy_split_model);

        auto const& strain_energy_tensile = _ip_data[ip].strain_energy_tensile;

        auto& ip_data = _ip_data[ip];
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
            (2 * N.transpose() * N * strain_energy_tensile) * w;

        local_rhs.noalias() -=
            (N.transpose() * N * d * 2 * strain_energy_tensile -
             local_pressure * dNdx.transpose() * N_u * u) *
            w;

        switch (_process_data.phasefield_model)
        {
            case PhaseFieldModel::AT1:
            {
                local_Jac.noalias() +=
                    gc * 0.75 * dNdx.transpose() * dNdx * ls * w;

                local_rhs.noalias() -=
                    gc *
                    (-0.375 * N.transpose() / ls +
                     0.75 * dNdx.transpose() * dNdx * ls * d) *
                    w;
                break;
            }
            case PhaseFieldModel::AT2:
            {
                local_Jac.noalias() +=
                    gc *
                    (N.transpose() * N / ls + dNdx.transpose() * dNdx * ls) * w;

                local_rhs.noalias() -=
                    gc *
                    ((N.transpose() * N / ls + dNdx.transpose() * dNdx * ls) *
                         d -
                     N.transpose() / ls) *
                    w;
                break;
            }
        }
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
                   [&](NumLib::LocalToGlobalIndexMap const& dof_table)
                   { return NumLib::getIndices(mesh_item_id, dof_table); });

    auto local_coupled_xs =
        getCoupledLocalSolutions(cpl_xs->coupled_xs, indices_of_processes);
    assert(local_coupled_xs.size() == displacement_size + phasefield_size);

    auto const d = Eigen::Map<PhaseFieldVector const>(
        &local_coupled_xs[phasefield_index], phasefield_size);
    auto const u = Eigen::Map<DeformationVector const>(
        &local_coupled_xs[displacement_index], displacement_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

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
                   [&](NumLib::LocalToGlobalIndexMap const& dof_table)
                   { return NumLib::getIndices(mesh_item_id, dof_table); });

    auto const local_coupled_xs =
        getCoupledLocalSolutions(cpl_xs->coupled_xs, indices_of_processes);
    assert(local_coupled_xs.size() == displacement_size + phasefield_size);

    auto const d = Eigen::Map<PhaseFieldVector const>(
        &local_coupled_xs[phasefield_index], phasefield_size);
    auto const u = Eigen::Map<DeformationVector const>(
        &local_coupled_xs[displacement_index], displacement_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    double element_elastic_energy = 0.0;
    double element_surface_energy = 0.0;
    double element_pressure_work = 0.0;

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;
        double const d_ip = N.dot(d);
        auto pressure_ip = _process_data.pressure;
        double const gc = _process_data.crack_resistance(t, x_position)[0];
        double const ls = _process_data.crack_length_scale(t, x_position)[0];

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

        element_elastic_energy += _ip_data[ip].elastic_energy * w;

        switch (_process_data.phasefield_model)
        {
            case PhaseFieldModel::AT1:
            {
                element_surface_energy +=
                    gc * 0.375 *
                    ((1 - d_ip) / ls + (dNdx * d).dot((dNdx * d)) * ls) * w;

                break;
            }
            case PhaseFieldModel::AT2:
            {
                element_surface_energy += 0.5 * gc *
                                          ((1 - d_ip) * (1 - d_ip) / ls +
                                           (dNdx * d).dot((dNdx * d)) * ls) *
                                          w;
                break;
            }
        }

        if (_process_data.crack_pressure)
        {
            element_pressure_work += pressure_ip * (N_u * u).dot(dNdx * d) * w;
        }
    }

#ifdef USE_PETSC
    int const n_all_nodes = indices_of_processes[1].size();
    int const n_regular_nodes = std::count_if(
        begin(indices_of_processes[1]), end(indices_of_processes[1]),
        [](GlobalIndexType const& index) { return index >= 0; });
    if (n_all_nodes != n_regular_nodes)
    {
        element_elastic_energy *=
            static_cast<double>(n_regular_nodes) / n_all_nodes;
        element_surface_energy *=
            static_cast<double>(n_regular_nodes) / n_all_nodes;
        element_pressure_work *=
            static_cast<double>(n_regular_nodes) / n_all_nodes;
    }
#endif  // USE_PETSC
    elastic_energy += element_elastic_energy;
    surface_energy += element_surface_energy;
    pressure_work += element_pressure_work;
}
}  // namespace PhaseField
}  // namespace ProcessLib
