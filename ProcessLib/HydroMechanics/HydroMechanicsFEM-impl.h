/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   HydroMechanicsFEM-impl.h
 *  Created on November 29, 2017, 2:03 PM
 */

#pragma once

#include "HydroMechanicsFEM.h"

#include "MaterialLib/SolidModels/KelvinVector.h"
#include "NumLib/Function/Interpolation.h"

namespace ProcessLib
{
namespace HydroMechanics
{
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void HydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                  ShapeFunctionPressure, IntegrationMethod,
                                  DisplacementDim>::
    assembleWithJacobianForPressureEquations(
        const double t, const std::vector<double>& local_xdot,
        const double /*dxdot_dx*/, const double /*dx_dx*/,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        const LocalCoupledSolutions& local_coupled_solutions)
{
    auto const& local_p = local_coupled_solutions.local_coupled_xs[0];
    auto const& local_u = local_coupled_solutions.local_coupled_xs[1];
    assert(local_p.size() == pressure_size);
    assert(local_u.size() == displacement_size);

    auto const local_matrix_size = local_p.size();
    auto local_rhs =
        MathLib::createZeroedVector<typename ShapeMatricesTypeDisplacement::
                                        template VectorType<pressure_size>>(
            local_b_data, local_matrix_size);

    SpatialPosition pos;
    pos.setElementID(this->_element.getID());

    auto p = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
        pressure_size> const>(local_p.data(), pressure_size);
    auto u =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_u.data(), displacement_size);

    auto p_dot =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            pressure_size> const>(local_xdot.data(), pressure_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            pressure_size, pressure_size>>(local_Jac_data, pressure_size,
                                           pressure_size);
    typename ShapeMatricesTypePressure::NodalMatrixType mass =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType laplace =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    double const& dt = _process_data.dt;

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;

        auto const& N_u = _ip_data[ip].N_u;
        auto const& dNdx_u = _ip_data[ip].dNdx_u;

        auto const& N_p = _ip_data[ip].N_p;
        auto const& dNdx_p = _ip_data[ip].dNdx_p;

        double const S = _process_data.specific_storage(t, x_position)[0];
        double const K_over_mu =
            _process_data.intrinsic_permeability(t, x_position)[0] /
            _process_data.fluid_viscosity(t, x_position)[0];
        auto const alpha_b = _process_data.biot_coefficient(t, x_position)[0];
        auto const rho_fr = _process_data.fluid_density(t, x_position)[0];

        laplace.noalias() += dNdx_p.transpose() * K_over_mu * dNdx_p * w;

        mass.noalias() += N_p.transpose() * S * N_p * w;

        auto const& b = _process_data.specific_body_force;
        local_rhs.noalias() += dNdx_p.transpose() * rho_fr * K_over_mu * b * w;

        auto& eps = _ip_data[ip].eps;
        auto const x_coord =
            interpolateXCoordinate<ShapeFunctionDisplacement,
                                   ShapeMatricesTypeDisplacement>(_element,
                                                                  N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        eps.noalias() = B * u;
        auto& eps_prev = _ip_data[ip].eps_prev;
        const double dv_dt =
            (Invariants::trace(eps) - Invariants::trace(eps_prev)) / dt;
        local_rhs.noalias() -= alpha_b * dv_dt * N_p * w;
    }
    local_Jac.noalias() = laplace + mass / dt;

    local_rhs.noalias() -= laplace * p + mass * p_dot;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void HydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                  ShapeFunctionPressure, IntegrationMethod,
                                  DisplacementDim>::
    assembleWithJacobianForDeformationEquations(
        const double t, const std::vector<double>& /*local_xdot*/,
        const double /*dxdot_dx*/, const double /*dx_dx*/,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        const LocalCoupledSolutions& local_coupled_solutions)
{
    auto const& local_p = local_coupled_solutions.local_coupled_xs[0];
    auto const& local_u = local_coupled_solutions.local_coupled_xs[1];
    assert(local_p.size() == pressure_size);
    assert(local_u.size() == displacement_size);

    auto u =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_u.data(), displacement_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size, displacement_size>>(
        local_Jac_data, displacement_size, displacement_size);

    auto local_rhs =
        MathLib::createZeroedVector<typename ShapeMatricesTypeDisplacement::
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

        auto const& N_u_op = _ip_data[ip].N_u_op;

        auto const& N_u = _ip_data[ip].N_u;
        auto const& dNdx_u = _ip_data[ip].dNdx_u;

        auto const& N_p = _ip_data[ip].N_p;

        auto const x_coord =
            interpolateXCoordinate<ShapeFunctionDisplacement,
                                   ShapeMatricesTypeDisplacement>(_element,
                                                                  N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        auto& eps = _ip_data[ip].eps;
        auto const& sigma_eff = _ip_data[ip].sigma_eff;

        auto const alpha = _process_data.biot_coefficient(t, x_position)[0];
        auto const rho_sr = _process_data.solid_density(t, x_position)[0];
        auto const rho_fr = _process_data.fluid_density(t, x_position)[0];
        auto const porosity = _process_data.porosity(t, x_position)[0];
        auto const& b = _process_data.specific_body_force;
        auto const& identity2 = MaterialLib::SolidModels::Invariants<
            KelvinVectorDimensions<DisplacementDim>::value>::identity2;

        eps.noalias() = B * u;

        auto C = _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u);

        local_Jac.noalias() += B.transpose() * C * B * w;

        double p_at_xi = 0.;
        NumLib::shapeFunctionInterpolate(local_p, N_p, p_at_xi);

        double const rho = rho_sr * (1 - porosity) + porosity * rho_fr;
        local_rhs.noalias() -=
            (B.transpose() * (sigma_eff - alpha * identity2 * p_at_xi) -
             N_u_op.transpose() * rho * b) *
            w;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void HydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                  ShapeFunctionPressure, IntegrationMethod,
                                  DisplacementDim>::
    assembleWithJacobianAndCoupling(
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
    if (local_coupled_solutions.process_id == 0)
    {
        assembleWithJacobianForPressureEquations(
            t, local_xdot, dxdot_dx, dx_dx, local_M_data, local_K_data,
            local_b_data, local_Jac_data, local_coupled_solutions);
        return;
    }

    // For the equations with deformation
    assembleWithJacobianForDeformationEquations(
        t, local_xdot, dxdot_dx, dx_dx, local_M_data, local_K_data,
        local_b_data, local_Jac_data, local_coupled_solutions);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void HydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                  ShapeFunctionPressure, IntegrationMethod,
                                  DisplacementDim>::
    postNonLinearSolverConcrete(std::vector<double> const& local_x,
                                double const t,
                                bool const use_monolithic_scheme)
{
    const int offset = use_monolithic_scheme ? displacement_index : 0;

    auto u =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x.data() + offset,
                                      displacement_size);
    double const& dt = _process_data.dt;
    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& N_u = _ip_data[ip].N_u;
        auto const& dNdx_u = _ip_data[ip].dNdx_u;

        auto const x_coord =
            interpolateXCoordinate<ShapeFunctionDisplacement,
                                   ShapeMatricesTypeDisplacement>(_element,
                                                                  N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        auto& eps = _ip_data[ip].eps;
        eps.noalias() = B * u;

        _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u);
    }
}

}  // namespace HydroMechanics
}  // namespace ProcessLib
