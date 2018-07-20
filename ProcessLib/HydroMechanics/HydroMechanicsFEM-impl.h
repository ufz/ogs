/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   HydroMechanicsFEM-impl.h
 *  Created on November 29, 2017, 2:03 PM
 */

#pragma once

#include "HydroMechanicsFEM.h"

#include "MathLib/KelvinVector.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
namespace HydroMechanics
{
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
HydroMechanicsLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                             IntegrationMethod, DisplacementDim>::
    HydroMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        HydroMechanicsProcessData<DisplacementDim>& process_data)
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
        initShapeMatrices<ShapeFunctionDisplacement,
                          ShapeMatricesTypeDisplacement, IntegrationMethod,
                          DisplacementDim>(e, is_axially_symmetric,
                                           _integration_method);

    auto const shape_matrices_p =
        initShapeMatrices<ShapeFunctionPressure, ShapeMatricesTypePressure,
                          IntegrationMethod, DisplacementDim>(
            e, is_axially_symmetric, _integration_method);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        // displacement (subscript u)
        _ip_data.emplace_back(*_process_data.material);
        auto& ip_data = _ip_data[ip];
        auto const& sm_u = shape_matrices_u[ip];
        _ip_data[ip].integration_weight =
            _integration_method.getWeightedPoint(ip).getWeight() *
            sm_u.integralMeasure * sm_u.detJ;

        // Initialize current time step values
        static const int kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
        ip_data.sigma_eff.setZero(kelvin_vector_size);
        ip_data.eps.setZero(kelvin_vector_size);

        // Previous time step values are not initialized and are set later.
        ip_data.eps_prev.resize(kelvin_vector_size);
        ip_data.sigma_eff_prev.resize(kelvin_vector_size);

        ip_data.N_u_op = ShapeMatricesTypeDisplacement::template MatrixType<
            DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                      displacement_size);
        for (int i = 0; i < DisplacementDim; ++i)
            ip_data.N_u_op
                .template block<1, displacement_size / DisplacementDim>(
                    i, i * displacement_size / DisplacementDim)
                .noalias() = sm_u.N;

        ip_data.N_u = sm_u.N;
        ip_data.dNdx_u = sm_u.dNdx;

        ip_data.N_p = shape_matrices_p[ip].N;
        ip_data.dNdx_p = shape_matrices_p[ip].dNdx;

        _secondary_data.N_u[ip] = shape_matrices_u[ip].N;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void HydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                  ShapeFunctionPressure, IntegrationMethod,
                                  DisplacementDim>::
    assembleWithJacobian(double const t, std::vector<double> const& local_x,
                         std::vector<double> const& local_xdot,
                         const double /*dxdot_dx*/, const double /*dx_dx*/,
                         std::vector<double>& /*local_M_data*/,
                         std::vector<double>& /*local_K_data*/,
                         std::vector<double>& local_rhs_data,
                         std::vector<double>& local_Jac_data)
{
    assert(local_x.size() == pressure_size + displacement_size);

    auto p = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
        pressure_size> const>(local_x.data() + pressure_index, pressure_size);

    auto u =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x.data() + displacement_index,
                                      displacement_size);

    auto p_dot =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            pressure_size> const>(local_xdot.data() + pressure_index,
                                  pressure_size);
    auto u_dot =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_xdot.data() + displacement_index,
                                      displacement_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size + pressure_size,
            displacement_size + pressure_size>>(
        local_Jac_data, displacement_size + pressure_size,
        displacement_size + pressure_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size + pressure_size>>(
        local_rhs_data, displacement_size + pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType laplace_p =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage_p =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        displacement_size, pressure_size>
        Kup = ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size, pressure_size>::Zero(displacement_size,
                                                    pressure_size);

    double const& dt = _process_data.dt;

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;

        auto const& N_u_op = _ip_data[ip].N_u_op;

        auto const& N_u = _ip_data[ip].N_u;
        auto const& dNdx_u = _ip_data[ip].dNdx_u;

        auto const& N_p = _ip_data[ip].N_p;
        auto const& dNdx_p = _ip_data[ip].dNdx_p;

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

        double const S = _process_data.specific_storage(t, x_position)[0];
        double const K_over_mu =
            _process_data.intrinsic_permeability(t, x_position)[0] /
            _process_data.fluid_viscosity(t, x_position)[0];
        auto const alpha = _process_data.biot_coefficient(t, x_position)[0];
        auto const rho_sr = _process_data.solid_density(t, x_position)[0];
        auto const rho_fr = _process_data.fluid_density(t, x_position)[0];
        auto const porosity = _process_data.porosity(t, x_position)[0];
        auto const& b = _process_data.specific_body_force;
        auto const& identity2 = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value>::identity2;

        //
        // displacement equation, displacement part
        //
        eps.noalias() = B * u;

        auto C = _ip_data[ip].updateConstitutiveRelation(
            t, x_position, dt, u, _process_data.reference_temperature);

        local_Jac
            .template block<displacement_size, displacement_size>(
                displacement_index, displacement_index)
            .noalias() += B.transpose() * C * B * w;

        double const rho = rho_sr * (1 - porosity) + porosity * rho_fr;
        local_rhs.template segment<displacement_size>(displacement_index)
            .noalias() -=
            (B.transpose() * sigma_eff - N_u_op.transpose() * rho * b) * w;

        //
        // displacement equation, pressure part
        //
        Kup.noalias() += B.transpose() * alpha * identity2 * N_p * w;

        //
        // pressure equation, pressure part.
        //
        laplace_p.noalias() += dNdx_p.transpose() * K_over_mu * dNdx_p * w;

        storage_p.noalias() += N_p.transpose() * S * N_p * w;

        local_rhs.template segment<pressure_size>(pressure_index).noalias() +=
            dNdx_p.transpose() * rho_fr * K_over_mu * b * w;

        //
        // pressure equation, displacement part.
        //
        // Reusing Kup.transpose().
    }
    // displacement equation, pressure part
    local_Jac
        .template block<displacement_size, pressure_size>(displacement_index,
                                                          pressure_index)
        .noalias() = -Kup;

    // pressure equation, pressure part.
    local_Jac
        .template block<pressure_size, pressure_size>(pressure_index,
                                                      pressure_index)
        .noalias() = laplace_p + storage_p / dt;

    // pressure equation, displacement part.
    local_Jac
        .template block<pressure_size, displacement_size>(pressure_index,
                                                          displacement_index)
        .noalias() = Kup.transpose() / dt;

    // pressure equation
    local_rhs.template segment<pressure_size>(pressure_index).noalias() -=
        laplace_p * p + storage_p * p_dot + Kup.transpose() * u_dot;

    // displacement equation
    local_rhs.template segment<displacement_size>(displacement_index)
        .noalias() += Kup * p;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const&
HydroMechanicsLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                             IntegrationMethod, DisplacementDim>::
    getIntPtDarcyVelocity(const double t,
                          GlobalVector const& current_solution,
                          NumLib::LocalToGlobalIndexMap const& dof_table,
                          std::vector<double>& cache) const
{
    auto const num_intpts = _ip_data.size();

    auto const indices = NumLib::getIndices(_element.getID(), dof_table);
    assert(!indices.empty());
    auto const local_x = current_solution.get(indices);

    cache.clear();
    auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, DisplacementDim, num_intpts);

    SpatialPosition pos;
    pos.setElementID(_element.getID());

    auto p = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
        pressure_size> const>(local_x.data() + pressure_index, pressure_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        double const K_over_mu =
            _process_data.intrinsic_permeability(t, x_position)[0] /
            _process_data.fluid_viscosity(t, x_position)[0];

        auto const rho_fr = _process_data.fluid_density(t, x_position)[0];
        auto const& b = _process_data.specific_body_force;

        // Compute the velocity
        auto const& dNdx_p = _ip_data[ip].dNdx_p;
        cache_matrix.col(ip).noalias() =
            -K_over_mu * dNdx_p * p - K_over_mu * rho_fr * b;
    }

    return cache;
}

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
        auto const& identity2 = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value>::identity2;

        eps.noalias() = B * u;

        auto C = _ip_data[ip].updateConstitutiveRelation(
            t, x_position, dt, u, _process_data.reference_temperature);

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
    const int displacement_offset =
        use_monolithic_scheme ? displacement_index : 0;

    auto u =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x.data() + displacement_offset,
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

        _ip_data[ip].updateConstitutiveRelation(
            t, x_position, dt, u, _process_data.reference_temperature);
    }
}

}  // namespace HydroMechanics
}  // namespace ProcessLib
