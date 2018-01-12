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

#include "PhaseFieldFEM.h"

namespace ProcessLib
{
namespace PhaseField
{
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
    if (local_coupled_solutions.process_id == 0)
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
        double const t, std::vector<double> const& local_xdot,
        const double dxdot_dx, const double dx_dx,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions) const
{
    auto const& local_d = local_coupled_solutions.local_coupled_xs[0];
    auto const& local_u = local_coupled_solutions.local_coupled_xs[1];
    assert(local_d.size() == phasefield_size);
    assert(local_u.size() == displacement_size);

    auto const local_matrix_size = local_u.size();
    auto d = Eigen::Map<
        typename ShapeMatricesType::template VectorType<phasefield_size> const>(
        local_d.data(), phasefield_size);

    auto u = Eigen::Map<typename ShapeMatricesType::template VectorType<
        displacement_size> const>(local_u.data(), displacement_size);

    // May not needed: auto d_dot = Eigen::Map<
    //    typename ShapeMatricesType::template VectorType<phasefield_size>
    //    const>(
    //    local_xdot.data(), phasefield_size);

    auto local_Jac = MathLib::createZeroedMatrix<JacobianMatrix>(
        local_Jac_data, local_matrix_size, local_matrix_size);

    auto local_rhs =
        MathLib::createZeroedVector<RhsVector>(local_b_data, local_matrix_size);

    typename ShapeMatricesType::template MatrixType<displacement_size,
                                                    phasefield_size>
        Kud;
    Kud.setZero(displacement_size, phasefield_size);

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
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void PhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                              DisplacementDim>::
    assembleWithJacobianPhaseFiledEquations(
        double const t, std::vector<double> const& local_xdot,
        const double dxdot_dx, const double dx_dx,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions) const
{
    auto const& local_d = local_coupled_solutions.local_coupled_xs[0];
    auto const& local_u = local_coupled_solutions.local_coupled_xs[1];
    assert(local_d.size() == phasefield_size);
    assert(local_u.size() == displacement_size);

    auto const local_matrix_size = local_d.size();
    auto d = Eigen::Map<
        typename ShapeMatricesType::template VectorType<phasefield_size> const>(
        local_d.data(), phasefield_size);

    auto u = Eigen::Map<typename ShapeMatricesType::template VectorType<
        displacement_size> const>(local_u.data(), displacement_size);

    auto d_dot = Eigen::Map<
        typename ShapeMatricesType::template VectorType<phasefield_size> const>(
        local_xdot.data(), phasefield_size);

    auto local_Jac = MathLib::createZeroedMatrix<JacobianMatrix>(
        local_Jac_data, local_matrix_size, local_matrix_size);

    auto local_rhs =
        MathLib::createZeroedVector<RhsVector>(local_b_data, local_matrix_size);

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

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;
    }
}

}  // namespace PhaseField
}  // namespace ProcessLib