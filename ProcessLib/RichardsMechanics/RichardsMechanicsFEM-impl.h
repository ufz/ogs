/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   RichardsMechanicsFEM-impl.h
 *  Created on November 29, 2017, 2:03 PM
 */

#pragma once

#include "RichardsMechanicsFEM.h"

#include "MathLib/KelvinVector.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
namespace RichardsMechanics
{
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                ShapeFunctionPressure, IntegrationMethod,
                                DisplacementDim>::
    RichardsMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        RichardsMechanicsProcessData<DisplacementDim>& process_data)
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
        _ip_data.emplace_back(*_process_data.solid_material);
        auto& ip_data = _ip_data[ip];
        auto const& sm_u = shape_matrices_u[ip];
        _ip_data[ip].integration_weight =
            _integration_method.getWeightedPoint(ip).getWeight() *
            sm_u.integralMeasure * sm_u.detJ;

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
void RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::assemble(double const t,
                               std::vector<double> const& local_x,
                               std::vector<double>& local_M_data,
                               std::vector<double>& local_K_data,
                               std::vector<double>& local_rhs_data)
{
    assert(local_x.size() == pressure_size + displacement_size);

    auto p_L =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            pressure_size> const>(local_x.data() + pressure_index,
                                  pressure_size);

    auto u =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x.data() + displacement_index,
                                      displacement_size);

    auto K = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size + pressure_size,
            displacement_size + pressure_size>>(
        local_K_data, displacement_size + pressure_size,
        displacement_size + pressure_size);

    auto M = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size + pressure_size,
            displacement_size + pressure_size>>(
        local_M_data, displacement_size + pressure_size,
        displacement_size + pressure_size);

    auto rhs = MathLib::createZeroedVector<
        typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size + pressure_size>>(
        local_rhs_data, displacement_size + pressure_size);

    double const& dt = _process_data.dt;
    auto const material_id =
        _process_data.flow_material->getMaterialID(_element.getID());

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

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N_p, p_cap_ip);

        auto& eps = _ip_data[ip].eps;
        auto& S_L = _ip_data[ip].saturation;

        auto const alpha = _process_data.biot_coefficient(t, x_position)[0];
        auto const rho_SR = _process_data.solid_density(t, x_position)[0];
        auto const K_SR = _process_data.solid_bulk_modulus(t, x_position)[0];
        auto const K_LR = _process_data.fluid_bulk_modulus(t, x_position)[0];
        auto const temperature = _process_data.temperature(t, x_position)[0];
        auto const porosity = _process_data.flow_material->getPorosity(
            material_id, t, x_position, -p_cap_ip, temperature, p_cap_ip);
        auto const rho_LR = _process_data.flow_material->getFluidDensity(
            -p_cap_ip, temperature);
        auto const& b = _process_data.specific_body_force;
        auto const& identity2 = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value>::identity2;

        S_L = _process_data.flow_material->getSaturation(
            material_id, t, x_position, -p_cap_ip, temperature, p_cap_ip);

        double const dS_L_dp_cap =
            _process_data.flow_material->getSaturationDerivative(
                material_id, t, x_position, -p_cap_ip, temperature, S_L);

        double const k_rel =
            _process_data.flow_material->getRelativePermeability(
                t, x_position, -p_cap_ip, temperature, S_L);
        auto const mu = _process_data.flow_material->getFluidViscosity(
            -p_cap_ip, temperature);

        double const K_intrinsic =
            _process_data.intrinsic_permeability(t, x_position)[0];
        double const K_over_mu = K_intrinsic / mu;

        //
        // displacement equation, displacement part
        //
        eps.noalias() = B * u;

        auto C = _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u,
                                                         temperature);

        //
        // displacement equation, displacement part
        //
        K.template block<displacement_size, displacement_size>(
             displacement_index, displacement_index)
            .noalias() += B.transpose() * C * B * w;

        double const rho = rho_SR * (1 - porosity) + porosity * rho_LR;
        rhs.template segment<displacement_size>(displacement_index).noalias() +=
            N_u_op.transpose() * rho * b * w;

        //
        // pressure equation, pressure part.
        //
        double const a0 = S_L * (alpha - porosity) / K_SR;
        // Volumetric average specific storage of the solid and fluid phases.
        double const specific_storage =
            dS_L_dp_cap * (p_cap_ip * a0 - porosity) +
            S_L * (porosity / K_LR + a0);
        M.template block<pressure_size, pressure_size>(pressure_index,
                                                       pressure_index)
            .noalias() += N_p.transpose() * rho_LR * specific_storage * N_p * w;

        K.template block<pressure_size, pressure_size>(pressure_index,
                                                       pressure_index)
            .noalias() +=
            dNdx_p.transpose() * rho_LR * k_rel * K_over_mu * dNdx_p * w;

        rhs.template segment<pressure_size>(pressure_index).noalias() +=
            dNdx_p.transpose() * rho_LR * rho_LR * k_rel * K_over_mu * b * w;

        //
        // displacement equation, pressure part
        //
        K.template block<displacement_size, pressure_size>(displacement_index,
                                                           pressure_index)
            .noalias() -= B.transpose() * alpha * S_L * identity2 * N_p * w;

        //
        // pressure equation, displacement part.
        //
        M.template block<pressure_size, displacement_size>(pressure_index,
                                                           displacement_index)
            .noalias() += N_p.transpose() * S_L * rho_LR * alpha *
                          identity2.transpose() * B * w;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
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

    auto p_L =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            pressure_size> const>(local_x.data() + pressure_index,
                                  pressure_size);

    auto u =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x.data() + displacement_index,
                                      displacement_size);

    auto p_L_dot =
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

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        pressure_size, displacement_size>
        Kpu = ShapeMatricesTypeDisplacement::template MatrixType<
            pressure_size, displacement_size>::Zero(pressure_size,
                                                    displacement_size);

    double const& dt = _process_data.dt;
    auto const material_id =
        _process_data.flow_material->getMaterialID(_element.getID());

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

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N_p, p_cap_ip);

        double p_cap_dot_ip;
        NumLib::shapeFunctionInterpolate(-p_L_dot, N_p, p_cap_dot_ip);

        typename ShapeMatricesTypeDisplacement::GlobalDimVectorType u_ip(
            DisplacementDim);
        for (int i = 0; i < u_ip.size(); ++i)
        {
            NumLib::shapeFunctionInterpolate(
                u.segment(i * ShapeFunctionDisplacement::NPOINTS,
                          ShapeFunctionDisplacement::NPOINTS),
                N_u, u_ip.coeffRef(i));
        }

        auto& eps = _ip_data[ip].eps;
        auto& S_L = _ip_data[ip].saturation;
        auto const& sigma_eff = _ip_data[ip].sigma_eff;

        auto const alpha = _process_data.biot_coefficient(t, x_position)[0];
        auto const rho_SR = _process_data.solid_density(t, x_position)[0];
        auto const K_SR = _process_data.solid_bulk_modulus(t, x_position)[0];
        auto const K_LR = _process_data.fluid_bulk_modulus(t, x_position)[0];
        auto const temperature = _process_data.temperature(t, x_position)[0];

        auto const porosity = _process_data.flow_material->getPorosity(
            material_id, t, x_position, -p_cap_ip, temperature, p_cap_ip);
        auto const rho_LR = _process_data.flow_material->getFluidDensity(
            -p_cap_ip, temperature);
        auto const& b = _process_data.specific_body_force;
        auto const& identity2 = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value>::identity2;

        S_L = _process_data.flow_material->getSaturation(
            material_id, t, x_position, -p_cap_ip, temperature, p_cap_ip);

        double const dS_L_dp_cap =
            _process_data.flow_material->getSaturationDerivative(
                material_id, t, x_position, -p_cap_ip, temperature, S_L);

        double const d2S_L_dp_cap_2 =
            _process_data.flow_material->getSaturationDerivative2(
                material_id, t, x_position, -p_cap_ip, temperature, S_L);

        double const k_rel =
            _process_data.flow_material->getRelativePermeability(
                t, x_position, -p_cap_ip, temperature, S_L);
        auto const mu = _process_data.flow_material->getFluidViscosity(
            -p_cap_ip, temperature);

        double const K_intrinsic =
            _process_data.intrinsic_permeability(t, x_position)[0];
        double const K_over_mu = K_intrinsic / mu;
        //
        // displacement equation, displacement part
        //
        eps.noalias() = B * u;

        auto C = _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u,
                                                         temperature);

        local_Jac
            .template block<displacement_size, displacement_size>(
                displacement_index, displacement_index)
            .noalias() += B.transpose() * C * B * w;

        double const rho = rho_SR * (1 - porosity) + porosity * rho_LR;
        local_rhs.template segment<displacement_size>(displacement_index)
            .noalias() -=
            (B.transpose() * sigma_eff - N_u_op.transpose() * rho * b) * w;

        //
        // displacement equation, pressure part
        //
        Kup.noalias() += B.transpose() * alpha * S_L * identity2 * N_p * w;

        /* For future implementation including swelling.
        double const dsigma_eff_dp_cap = -K_intrinsic * m_swell * n *
                                         std::pow(S_L, n - 1) * dS_L_dp_cap *
                                         identity2;
        local_Jac
            .template block<displacement_size, pressure_size>(
                displacement_index, pressure_index)
            .noalias() -= B.transpose() * dsigma_eff_dp_cap * N_p * w;
        */

        local_Jac
            .template block<displacement_size, pressure_size>(
                displacement_index, pressure_index)
            .noalias() -= B.transpose() * alpha *
                          (S_L + p_cap_ip * dS_L_dp_cap) * identity2 * N_p * w;

        local_Jac
            .template block<displacement_size, pressure_size>(
                displacement_index, pressure_index)
            .noalias() +=
            N_u_op.transpose() * porosity * rho_LR * dS_L_dp_cap * b * N_p * w;
        //
        // pressure equation, displacement part.
        //
        Kpu.noalias() += N_p.transpose() * S_L * rho_LR * alpha *
                         identity2.transpose() * B * w;

        //
        // pressure equation, pressure part.
        //
        laplace_p.noalias() +=
            dNdx_p.transpose() * rho_LR * k_rel * K_over_mu * dNdx_p * w;

        double const a0 = (alpha - porosity) / K_SR;
        double const specific_storage =
            dS_L_dp_cap * (p_cap_ip * S_L * a0 - porosity) +
            S_L * (porosity / K_LR + S_L * a0);

        double const dspecific_storage_dp_cap =
            d2S_L_dp_cap_2 * (p_cap_ip * S_L * a0 - porosity) +
            dS_L_dp_cap *
                (porosity / K_LR + a0 * 3 * S_L + dS_L_dp_cap * p_cap_ip * a0);

        storage_p.noalias() +=
            N_p.transpose() * rho_LR * specific_storage * N_p * w;

        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() += N_p.transpose() * rho_LR * p_cap_dot_ip *
                          dspecific_storage_dp_cap * N_p * w;

        /* In the derivation there is a div(du/dt) term in the Jacobian, but
         * this implementation increases the total runtime by 1%. Maybe a very
         * large step is needed to see the increase of efficiency.
        double div_u_dot = 0;
        for (int i = 0; i < DisplacementDim; ++i)
        {
            div_u_dot +=
                (dNdx_u *
                 u_dot.template segment<ShapeFunctionDisplacement::NPOINTS>(
                     i * ShapeFunctionDisplacement::NPOINTS))[i];
        }
        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() -= N_p.transpose() * rho_LR * dS_L_dp_cap * alpha *
                          div_u_dot * N_p * w;
         */

        double const dk_rel_dS_l =
            _process_data.flow_material->getRelativePermeabilityDerivative(
                t, x_position, -p_cap_ip, temperature, S_L);
        typename ShapeMatricesTypeDisplacement::GlobalDimVectorType const
            grad_p_cap = -dNdx_p * p_L;
        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() += dNdx_p.transpose() * rho_LR * K_over_mu * grad_p_cap *
                          dk_rel_dS_l * dS_L_dp_cap * N_p * w;

        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() += dNdx_p.transpose() * rho_LR * rho_LR * K_over_mu * b *
                          dk_rel_dS_l * dS_L_dp_cap * N_p * w;

        local_rhs.template segment<pressure_size>(pressure_index).noalias() +=
            dNdx_p.transpose() * rho_LR * rho_LR * k_rel * K_over_mu * b * w;
    }

    // pressure equation, pressure part.
    local_Jac
        .template block<pressure_size, pressure_size>(pressure_index,
                                                      pressure_index)
        .noalias() += laplace_p + storage_p / dt;

    // pressure equation, displacement part.
    local_Jac
        .template block<pressure_size, displacement_size>(pressure_index,
                                                          displacement_index)
        .noalias() = Kpu / dt;

    // pressure equation
    local_rhs.template segment<pressure_size>(pressure_index).noalias() -=
        laplace_p * p_L + storage_p * p_L_dot + Kpu * u_dot;

    // displacement equation
    local_rhs.template segment<displacement_size>(displacement_index)
        .noalias() += Kup * p_L;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const& RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::
    getIntPtSigma(const double /*t*/,
                  GlobalVector const& /*current_solution*/,
                  NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                  std::vector<double>& cache) const
{
    static const int kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    auto const num_intpts = _ip_data.size();

    cache.clear();
    auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, kelvin_vector_size, num_intpts);

    for (unsigned ip = 0; ip < num_intpts; ++ip)
    {
        auto const& sigma = _ip_data[ip].sigma_eff;
        cache_mat.col(ip) =
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma);
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const& RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::
    getIntPtEpsilon(const double /*t*/,
                    GlobalVector const& /*current_solution*/,
                    NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                    std::vector<double>& cache) const
{
    auto const kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    auto const num_intpts = _ip_data.size();

    cache.clear();
    auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, kelvin_vector_size, num_intpts);

    for (unsigned ip = 0; ip < num_intpts; ++ip)
    {
        auto const& eps = _ip_data[ip].eps;
        cache_mat.col(ip) =
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(eps);
    }

    return cache;
}
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const& RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::getIntPtDarcyVelocity(const double t,
                                            GlobalVector const&
                                                current_solution,
                                            NumLib::LocalToGlobalIndexMap const&
                                                dof_table,
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

    auto p_L =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            pressure_size> const>(local_x.data() + pressure_index,
                                  pressure_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& N_p = _ip_data[ip].N_p;

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N_p, p_cap_ip);

        auto const temperature = _process_data.temperature(t, x_position)[0];
        auto const mu = _process_data.flow_material->getFluidViscosity(
            -p_cap_ip, temperature);
        double const K_over_mu =
            _process_data.intrinsic_permeability(t, x_position)[0] / mu;
        auto const rho_LR = _process_data.flow_material->getFluidDensity(
            -p_cap_ip, temperature);
        auto const& b = _process_data.specific_body_force;

        // Compute the velocity
        auto const& dNdx_p = _ip_data[ip].dNdx_p;
        cache_matrix.col(ip).noalias() =
            -K_over_mu * dNdx_p * p_L - K_over_mu * rho_LR * b;
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const& RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::
    getIntPtSaturation(const double /*t*/,
                       GlobalVector const& /*current_solution*/,
                       NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                       std::vector<double>& cache) const
{
    auto const num_intpts = _ip_data.size();

    cache.clear();
    auto cache_mat = MathLib::createZeroedMatrix<
        Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>>(cache, 1,
                                                                   num_intpts);

    for (unsigned ip = 0; ip < num_intpts; ++ip)
    {
        cache_mat[ip] = _ip_data[ip].saturation;
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, IntegrationMethod,
                                     DisplacementDim>::
    assembleWithJacobianForPressureEquations(
        const double /*t*/, const std::vector<double>& /*local_xdot*/,
        const double /*dxdot_dx*/, const double /*dx_dx*/,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& /*local_b_data*/,
        std::vector<double>& /*local_Jac_data*/,
        const LocalCoupledSolutions& /*local_coupled_solutions*/)
{
    OGS_FATAL("RichardsMechanics; The staggered scheme is not implemented.");
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, IntegrationMethod,
                                     DisplacementDim>::
    assembleWithJacobianForDeformationEquations(
        const double /*t*/, const std::vector<double>& /*local_xdot*/,
        const double /*dxdot_dx*/, const double /*dx_dx*/,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& /*local_b_data*/,
        std::vector<double>& /*local_Jac_data*/,
        const LocalCoupledSolutions& /*local_coupled_solutions*/)
{
    OGS_FATAL("RichardsMechanics; The staggered scheme is not implemented.");
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
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
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
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
        auto const temperature = _process_data.temperature(t, x_position)[0];

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

        _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u,
                                                temperature);
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, IntegrationMethod,
                                     DisplacementDim>::
    computeSecondaryVariableConcrete(double const t,
                                     std::vector<double> const& local_x)
{
    auto p_L =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            pressure_size> const>(local_x.data() + pressure_index,
                                  pressure_size);
    auto const material_id =
        _process_data.flow_material->getMaterialID(_element.getID());

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    double saturation_avg = 0;

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& N_p = _ip_data[ip].N_p;
        auto const temperature = _process_data.temperature(t, x_position)[0];

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N_p, p_cap_ip);
        auto& S_L = _ip_data[ip].saturation;
        S_L = _process_data.flow_material->getSaturation(
            material_id, t, x_position, -p_cap_ip, temperature, p_cap_ip);
        saturation_avg += S_L;
    }
    saturation_avg /= n_integration_points;

    (*_process_data.element_saturation)[_element.getID()] = saturation_avg;
}

}  // namespace RichardsMechanics
}  // namespace ProcessLib
