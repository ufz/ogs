/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   PhaseFieldFEM-impl.h
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
    assert(local_matrix_size == phasefield_size + displacement_size);

    auto d = Eigen::Map<
        typename ShapeMatricesType::template VectorType<phasefield_size> const>(
        local_x.data() + phasefield_index, phasefield_size);

    auto u = Eigen::Map<typename ShapeMatricesType::template VectorType<
        displacement_size> const>(local_x.data() + displacement_index,
                                  displacement_size);

    auto d_dot = Eigen::Map<
        typename ShapeMatricesType::template VectorType<phasefield_size> const>(
        local_xdot.data() + phasefield_index, phasefield_size);

    auto local_Jac = MathLib::createZeroedMatrix<JacobianMatrix>(
        local_Jac_data, local_matrix_size, local_matrix_size);

    auto local_rhs = MathLib::createZeroedVector<RhsVector>(local_rhs_data,
                                                            local_matrix_size);

    typename ShapeMatricesType::template MatrixType<displacement_size,
                                                    phasefield_size>
        Kud;
    Kud.setZero(displacement_size, phasefield_size);

    typename ShapeMatricesType::template MatrixType<phasefield_size,
                                                    displacement_size>
        Kdu;
    Kdu.setZero(phasefield_size, displacement_size);

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

        auto const& C_tensile = _ip_data[ip].C_tensile;
        auto const& C_compressive = _ip_data[ip].C_compressive;

        auto& eps = _ip_data[ip].eps;
        auto const& strain_energy_tensile = _ip_data[ip].strain_energy_tensile;
        auto const& sigma_tensile = _ip_data[ip].sigma_tensile;

        auto& history_variable = _ip_data[ip].history_variable;
        auto& history_variable_prev = _ip_data[ip].history_variable_prev;

        auto const& sigma_real = _ip_data[ip].sigma_real;

        // Kdd_1 defines one term which both used in Kdd and local_rhs for
        // phase field
        double const gc = _process_data.crack_resistance(t, x_position)[0];
        double const ls = _process_data.crack_length_scale(t, x_position)[0];
        typename ShapeMatricesType::NodalMatrixType const Kdd_1 =
            dNdx.transpose() * 2 * gc * ls * dNdx;

        //
        // displacement equation, displacement part
        //
        double const k = _process_data.residual_stiffness(t, x_position)[0];
        double const d_ip = N.dot(d);
        double const degradation = d_ip * d_ip * (1 - k) + k;
        eps.noalias() = B * u;
        _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u,
                                                degradation);

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

        auto const rho_sr = _process_data.solid_density(t, x_position)[0];
        auto const& b = _process_data.specific_body_force;
        local_rhs.template block<displacement_size, 1>(displacement_index, 0)
            .noalias() -=
            (B.transpose() * sigma_real - N_u.transpose() * rho_sr * b) * w;

        //
        // displacement equation, phasefield part
        //

        // Temporary storage used in the Kud and for potential reuse in the
        // Kdu matrix.
        auto const Kud_ip_contribution =
            (B.transpose() * 2 * d_ip * sigma_tensile * N * w).eval();

        Kud.noalias() += Kud_ip_contribution;

        if (history_variable_prev < strain_energy_tensile)
        {
            history_variable = strain_energy_tensile;
            Kdu.noalias() += Kud_ip_contribution.transpose();
        }
        else
        {
            history_variable = history_variable_prev;
        }

        //
        // phasefield equation, phasefield part.
        //
        Kdd.noalias() += (Kdd_1 + N.transpose() * 2 * history_variable * N +
                          N.transpose() * 0.5 * gc / ls * N) *
                         w;
        double const M = _process_data.kinetic_coefficient(t, x_position)[0];
        double const d_dot_ip = N.dot(d_dot);

        local_rhs.template segment<phasefield_size>(phasefield_index)
            .noalias() -= (N.transpose() * d_dot_ip / M + Kdd_1 * d +
                           N.transpose() * d_ip * 2 * history_variable -
                           N.transpose() * 0.5 * gc / ls * (1 - d_ip)) *
                          w;

        Ddd.noalias() += N.transpose() / M * N * w;
    }
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
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void PhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                              DisplacementDim>::
    assembleWithJacobianForStaggeredScheme(
        double const t, std::vector<double> const& local_xdot,
        const double dxdot_dx, const double dx_dx,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions)
{
    // For the equations with phase field.
    if (local_coupled_solutions.process_id == 1)
    {
        assembleWithJacobianPhaseFiledEquations(
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
void PhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                              DisplacementDim>::
    assembleWithJacobianForDeformationEquations(
        double const t, std::vector<double> const& /*local_xdot*/,
        const double /*dxdot_dx*/, const double /*dx_dx*/,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions)
{
    using DeformationVector =
        typename ShapeMatricesType::template VectorType<displacement_size>;
    using DeformationMatrix =
        typename ShapeMatricesType::template MatrixType<displacement_size,
                                                        displacement_size>;
    using PhaseFieldVector =
        typename ShapeMatricesType::template VectorType<phasefield_size>;

    auto const& local_u = local_coupled_solutions.local_coupled_xs[0];
    auto const& local_d = local_coupled_solutions.local_coupled_xs[1];
    assert(local_u.size() == displacement_size);
    assert(local_d.size() == phasefield_size);

    auto const local_matrix_size = local_u.size();
    auto d =
        Eigen::Map<PhaseFieldVector const>(local_d.data(), phasefield_size);

    auto u =
        Eigen::Map<DeformationVector const>(local_u.data(), displacement_size);

    auto local_Jac = MathLib::createZeroedMatrix<DeformationMatrix>(
        local_Jac_data, local_matrix_size, local_matrix_size);

    auto local_rhs = MathLib::createZeroedVector<DeformationVector>(
        local_b_data, local_matrix_size);

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());
    double const& dt = _process_data.dt;

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
        eps.noalias() = B * u;
        double const k = _process_data.residual_stiffness(t, x_position)[0];
        double const d_ip = N.dot(d);
        double const degradation = d_ip * d_ip * (1 - k) + k;
        _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u,
                                                degradation);

        auto& sigma_real = _ip_data[ip].sigma_real;
        auto const& C_tensile = _ip_data[ip].C_tensile;
        auto const& C_compressive = _ip_data[ip].C_compressive;

        typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                        displacement_size>
            N_u = ShapeMatricesType::template MatrixType<
                DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                          displacement_size);

        for (int i = 0; i < DisplacementDim; ++i)
            N_u.template block<1, displacement_size / DisplacementDim>(
                   i, i * displacement_size / DisplacementDim)
                .noalias() = N;

        auto const rho_sr = _process_data.solid_density(t, x_position)[0];
        auto const& b = _process_data.specific_body_force;

        auto pressure = _ip_data[ip].pressure;

        local_rhs.noalias() -=
            (B.transpose() * sigma_real - N_u.transpose() * rho_sr * b -
             pressure * N_u.transpose() * dNdx * d) *
            w;

        local_Jac.noalias() +=
            B.transpose() * (degradation * C_tensile + C_compressive) * B * w;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void PhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                              DisplacementDim>::
    assembleWithJacobianPhaseFiledEquations(
        double const t, std::vector<double> const& /*local_xdot*/,
        const double /*dxdot_dx*/, const double /*dx_dx*/,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions)
{
    using DeformationVector =
        typename ShapeMatricesType::template VectorType<displacement_size>;
    using PhaseFieldVector =
        typename ShapeMatricesType::template VectorType<phasefield_size>;
    using PhaseFieldMatrix =
        typename ShapeMatricesType::template MatrixType<phasefield_size,
                                                        phasefield_size>;

    auto const& local_u = local_coupled_solutions.local_coupled_xs[0];
    auto const& local_d = local_coupled_solutions.local_coupled_xs[1];
    assert(local_u.size() == displacement_size);
    assert(local_d.size() == phasefield_size);

    auto const local_matrix_size = local_d.size();
    auto d =
        Eigen::Map<PhaseFieldVector const>(local_d.data(), phasefield_size);
    auto u =
        Eigen::Map<DeformationVector const>(local_u.data(), displacement_size);

    auto local_Jac = MathLib::createZeroedMatrix<PhaseFieldMatrix>(
        local_Jac_data, local_matrix_size, local_matrix_size);
    auto local_rhs = MathLib::createZeroedVector<PhaseFieldVector>(
        local_b_data, local_matrix_size);

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        double const gc = _process_data.crack_resistance(t, x_position)[0];
        double const ls = _process_data.crack_length_scale(t, x_position)[0];
        auto const& strain_energy_tensile = _ip_data[ip].strain_energy_tensile;

        auto& ip_data = _ip_data[ip];
        ip_data.strain_energy_tensile = strain_energy_tensile;

        typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                        displacement_size>
            N_u = ShapeMatricesType::template MatrixType<
                DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                          displacement_size);

        for (int i = 0; i < DisplacementDim; ++i)
            N_u.template block<1, displacement_size / DisplacementDim>(
                   i, i * displacement_size / DisplacementDim)
                .noalias() = N;

        auto pressure = _ip_data[ip].pressure;

        local_Jac.noalias() +=
            (2 * N.transpose() * N * strain_energy_tensile +
             gc * (N.transpose() * N / ls + dNdx.transpose() * dNdx * ls)) *
            w;

        local_rhs.noalias() -=
            (N.transpose() * N * d * 2 * strain_energy_tensile +
             gc * ((N.transpose() * N / ls + dNdx.transpose() * dNdx * ls) * d -
                   N.transpose() / ls) -
             pressure * dNdx.transpose() * N_u * u) *
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
                         bool const /*use_monolithic_scheme*/,
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
        getCurrentLocalSolutions(*cpl_xs, indices_of_processes);
    assert(local_coupled_xs.size() == 2);

    auto const& local_u = local_coupled_xs[0];
    auto const& local_d = local_coupled_xs[1];

    assert(local_u.size() == displacement_size);
    assert(local_d.size() == phasefield_size);

    auto d = Eigen::Map<
        typename ShapeMatricesType::template VectorType<phasefield_size> const>(
        local_d.data(), phasefield_size);

    auto u = Eigen::Map<typename ShapeMatricesType::template VectorType<
        displacement_size> const>(local_u.data(), displacement_size);

    SpatialPosition x_position;
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
            N_u.template block<1, displacement_size / DisplacementDim>(
                   i, i * displacement_size / DisplacementDim)
                .noalias() = N;

        crack_volume += (N_u * u).dot(dNdx * d) * w;
    }
}
}  // namespace PhaseField
}  // namespace ProcessLib
