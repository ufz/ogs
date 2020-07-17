/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on November 29, 2017, 2:03 PM
 */

#pragma once

#include "HydroMechanicsFEM.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/MPL/Utils/GetSymmetricTensor.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"

namespace ProcessLib
{
namespace HydroMechanics
{
namespace MPL = MaterialPropertyLib;

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
        NumLib::initShapeMatrices<ShapeFunctionDisplacement,
                                  ShapeMatricesTypeDisplacement,
                                  DisplacementDim>(e, is_axially_symmetric,
                                                   _integration_method);

    auto const shape_matrices_p =
        NumLib::initShapeMatrices<ShapeFunctionPressure,
                                  ShapeMatricesTypePressure, DisplacementDim>(
            e, is_axially_symmetric, _integration_method);

    auto const& solid_material =
        MaterialLib::Solids::selectSolidConstitutiveRelation(
            _process_data.solid_materials, _process_data.material_ids,
            e.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        _ip_data.emplace_back(solid_material);
        auto& ip_data = _ip_data[ip];
        auto const& sm_u = shape_matrices_u[ip];
        ip_data.integration_weight =
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
        {
            ip_data.N_u_op
                .template block<1, displacement_size / DisplacementDim>(
                    i, i * displacement_size / DisplacementDim)
                .noalias() = sm_u.N;
        }

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
    assembleWithJacobian(double const t, double const dt,
                         std::vector<double> const& local_x,
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

    typename ShapeMatricesTypePressure::NodalMatrixType add_p_derivative =
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

    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material =
        *_process_data.solid_materials[0];

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    auto const& b = _process_data.specific_body_force;
    auto const& medium = _process_data.media_map->getMedium(_element.getID());
    auto const& solid = medium->phase("Solid");
    auto const& gas = medium->phase("Gas");
    MPL::VariableArray vars;

    auto const T_ref =
        medium->property(MPL::PropertyType::reference_temperature)
            .template value<double>(vars, x_position, t, dt);
    vars[static_cast<int>(MPL::Variable::temperature)] = T_ref;

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value>::
        identity2;

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
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                _element, N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        auto& eps = _ip_data[ip].eps;
        auto const& sigma_eff = _ip_data[ip].sigma_eff;

        double const p_int_pt = N_p.dot(p);
        vars[static_cast<int>(MPL::Variable::phase_pressure)] = p_int_pt;

        auto const K_S = solid_material.getBulkModulus(t, x_position);

        auto const alpha = solid.property(MPL::PropertyType::biot_coefficient)
                               .template value<double>(vars, x_position, t, dt);

        // For stress dependent permeability.
        vars[static_cast<int>(MPL::Variable::total_stress)].emplace<SymmetricTensor>(
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                (_ip_data[ip].sigma_eff - alpha * identity2 * p_int_pt)
                    .eval()));

        auto const K = MPL::formEigenTensor<DisplacementDim>(
            medium->property(MPL::PropertyType::permeability)
                .value(vars, x_position, t, dt));
        auto const rho_sr =
            solid.property(MPL::PropertyType::density)
                .template value<double>(vars, x_position, t, dt);
        auto const porosity =
            solid.property(MPL::PropertyType::porosity)
                .template value<double>(vars, x_position, t, dt);

        auto const mu = gas.property(MPL::PropertyType::viscosity)
                            .template value<double>(vars, x_position, t, dt);
        auto const rho_fr =
            gas.property(MPL::PropertyType::density)
                .template value<double>(vars, x_position, t, dt);
        auto const beta_p =
            gas.property(MPL::PropertyType::density)
                .template dValue<double>(vars, MPL::Variable::phase_pressure,
                                         x_position, t, dt) /
            rho_fr;

        auto const K_over_mu = K / mu;

        //
        // displacement equation, displacement part
        //
        eps.noalias() = B * u;
        vars[static_cast<int>(MPL::Variable::strain)]
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps);

        auto C = _ip_data[ip].updateConstitutiveRelation(vars, t, x_position,
                                                         dt, u, T_ref);

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
        laplace_p.noalias() +=
            rho_fr * dNdx_p.transpose() * K_over_mu * dNdx_p * w;

        storage_p.noalias() +=
            rho_fr * N_p.transpose() * N_p * w *
            ((alpha - porosity) * (1.0 - alpha) / K_S + porosity * beta_p);

        // density dependence on pressure evaluated for Darcy-term,
        // for laplace and storage terms this dependence is neglected
        add_p_derivative.noalias() += rho_fr * beta_p * dNdx_p.transpose() *
                                      K_over_mu *
                                      (dNdx_p * p - 2.0 * rho_fr * b) * N_p * w;

        local_rhs.template segment<pressure_size>(pressure_index).noalias() +=
            dNdx_p.transpose() * rho_fr * rho_fr * K_over_mu * b * w;

        //
        // pressure equation, displacement part.
        //
        Kpu.noalias() +=
            rho_fr * alpha * N_p.transpose() * identity2.transpose() * B * w;
    }
    // displacement equation, pressure part
    local_Jac
        .template block<displacement_size, pressure_size>(displacement_index,
                                                          pressure_index)
        .noalias() = -Kup;

    if (_process_data.mass_lumping)
    {
        storage_p = storage_p.colwise().sum().eval().asDiagonal();
    }

    // pressure equation, pressure part.
    local_Jac
        .template block<pressure_size, pressure_size>(pressure_index,
                                                      pressure_index)
        .noalias() = laplace_p + storage_p / dt + add_p_derivative;

    // pressure equation, displacement part.
    local_Jac
        .template block<pressure_size, displacement_size>(pressure_index,
                                                          displacement_index)
        .noalias() = Kpu / dt;

    // pressure equation
    local_rhs.template segment<pressure_size>(pressure_index).noalias() -=
        laplace_p * p + storage_p * p_dot + Kpu * u_dot;

    // displacement equation
    local_rhs.template segment<displacement_size>(displacement_index)
        .noalias() += Kup * p;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const&
HydroMechanicsLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                             IntegrationMethod, DisplacementDim>::
    getIntPtDarcyVelocity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const
{
    int const hydraulic_process_id = _process_data.hydraulic_process_id;
    auto const indices =
        NumLib::getIndices(_element.getID(), *dof_table[hydraulic_process_id]);
    assert(!indices.empty());
    auto const local_x = x[hydraulic_process_id]->get(indices);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    cache.clear();
    auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, DisplacementDim, n_integration_points);

    auto p = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
        pressure_size> const>(local_x.data() + pressure_index, pressure_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    auto const& medium = _process_data.media_map->getMedium(_element.getID());
    auto const& gas = medium->phase("Gas");
    auto const& solid = medium->phase("Solid");
    MPL::VariableArray vars;

    // TODO (naumov) Temporary value not used by current material models. Need
    // extension of secondary variables interface.
    double const dt = std::numeric_limits<double>::quiet_NaN();
    vars[static_cast<int>(MPL::Variable::temperature)] =
        medium->property(MPL::PropertyType::reference_temperature)
            .template value<double>(vars, x_position, t, dt);

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value>::
        identity2;

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        double const p_int_pt = _ip_data[ip].N_p.dot(p);
        vars[static_cast<int>(MPL::Variable::phase_pressure)] = p_int_pt;

        auto const alpha = solid.property(MPL::PropertyType::biot_coefficient)
                               .template value<double>(vars, x_position, t, dt);

        // For stress dependent permeability.
        vars[static_cast<int>(MPL::Variable::total_stress)].emplace<SymmetricTensor>(
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                (_ip_data[ip].sigma_eff - alpha * identity2 * p_int_pt)
                    .eval()));
        auto const K = MPL::formEigenTensor<DisplacementDim>(
            medium->property(MPL::PropertyType::permeability)
                .value(vars, x_position, t, dt));

        auto const mu = gas.property(MPL::PropertyType::viscosity)
                            .template value<double>(vars, x_position, t, dt);
        auto const rho_fr =
            gas.property(MPL::PropertyType::density)
                .template value<double>(vars, x_position, t, dt);

        auto const K_over_mu = K / mu;

        auto const& b = _process_data.specific_body_force;

        // Compute the velocity
        auto const& dNdx_p = _ip_data[ip].dNdx_p;
        cache_matrix.col(ip).noalias() =
            -K_over_mu * dNdx_p * p + K_over_mu * rho_fr * b;
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void HydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                  ShapeFunctionPressure, IntegrationMethod,
                                  DisplacementDim>::
    assembleWithJacobianForPressureEquations(
        const double t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_xdot, std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data)
{
    auto local_rhs =
        MathLib::createZeroedVector<typename ShapeMatricesTypeDisplacement::
                                        template VectorType<pressure_size>>(
            local_b_data, pressure_size);

    ParameterLib::SpatialPosition pos;
    pos.setElementID(this->_element.getID());

    auto const p = local_x.template segment<pressure_size>(pressure_index);
    auto const u =
        local_x.template segment<displacement_size>(displacement_index);

    auto const p_dot =
        local_xdot.template segment<pressure_size>(pressure_index);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            pressure_size, pressure_size>>(local_Jac_data, pressure_size,
                                           pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType laplace =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType add_p_derivative =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material =
        *_process_data.solid_materials[0];

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    auto const& medium = _process_data.media_map->getMedium(_element.getID());
    auto const& solid = medium->phase("Solid");
    auto const& gas = medium->phase("Gas");
    MPL::VariableArray vars;

    vars[static_cast<int>(MPL::Variable::temperature)] =
        medium->property(MPL::PropertyType::reference_temperature)
            .template value<double>(vars, x_position, t, dt);
    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value>::
        identity2;

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;

        auto const& N_u = _ip_data[ip].N_u;
        auto const& dNdx_u = _ip_data[ip].dNdx_u;

        auto const& N_p = _ip_data[ip].N_p;
        auto const& dNdx_p = _ip_data[ip].dNdx_p;

        double const p_int_pt = N_p.dot(p);
        vars[static_cast<int>(MPL::Variable::phase_pressure)] = p_int_pt;

        auto const K_S = solid_material.getBulkModulus(t, x_position);

        auto const alpha_b =
            solid.property(MPL::PropertyType::biot_coefficient)
                .template value<double>(vars, x_position, t, dt);

        // For stress dependent permeability.
        vars[static_cast<int>(MPL::Variable::total_stress)].emplace<SymmetricTensor>(
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                (_ip_data[ip].sigma_eff - alpha_b * identity2 * p_int_pt)
                    .eval()));

        auto const K = MPL::formEigenTensor<DisplacementDim>(
            medium->property(MPL::PropertyType::permeability)
                .value(vars, x_position, t, dt));
        auto const porosity =
            solid.property(MPL::PropertyType::porosity)
                .template value<double>(vars, x_position, t, dt);

        auto const mu = gas.property(MPL::PropertyType::viscosity)
                            .template value<double>(vars, x_position, t, dt);
        auto const rho_fr =
            gas.property(MPL::PropertyType::density)
                .template value<double>(vars, x_position, t, dt);
        auto const beta_p =
            gas.property(MPL::PropertyType::density)
                .template dValue<double>(vars, MPL::Variable::phase_pressure,
                                         x_position, t, dt) /
            rho_fr;

        auto const K_over_mu = K / mu;

        laplace.noalias() +=
            rho_fr * dNdx_p.transpose() * K_over_mu * dNdx_p * w;

        storage.noalias() +=
            rho_fr * N_p.transpose() * N_p * w *
            ((alpha_b - porosity) * (1.0 - alpha_b) / K_S + porosity * beta_p);

        auto const& b = _process_data.specific_body_force;
        local_rhs.noalias() +=
            dNdx_p.transpose() * rho_fr * rho_fr * K_over_mu * b * w;

        // density dependence on pressure evaluated for Darcy-term,
        // for laplace and storage terms this dependence is neglected (as is
        // done for monolithic too)
        add_p_derivative.noalias() += rho_fr * beta_p * dNdx_p.transpose() *
                                      K_over_mu *
                                      (dNdx_p * p - 2.0 * rho_fr * b) * N_p * w;

        auto& eps = _ip_data[ip].eps;
        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                _element, N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        eps.noalias() = B * u;
        auto& eps_prev = _ip_data[ip].eps_prev;
        const double dv_dt =
            (Invariants::trace(eps) - Invariants::trace(eps_prev)) / dt;
        local_rhs.noalias() -= rho_fr * alpha_b * dv_dt * N_p * w;
    }
    local_Jac.noalias() = laplace + storage / dt + add_p_derivative;

    local_rhs.noalias() -= laplace * p + storage * p_dot;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void HydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                  ShapeFunctionPressure, IntegrationMethod,
                                  DisplacementDim>::
    assembleWithJacobianForDeformationEquations(
        const double t, double const dt, Eigen::VectorXd const& local_x,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    auto const p = local_x.template segment<pressure_size>(pressure_index);
    auto const u =
        local_x.template segment<displacement_size>(displacement_index);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size, displacement_size>>(
        local_Jac_data, displacement_size, displacement_size);

    auto local_rhs =
        MathLib::createZeroedVector<typename ShapeMatricesTypeDisplacement::
                                        template VectorType<displacement_size>>(
            local_b_data, displacement_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    auto const& medium = _process_data.media_map->getMedium(_element.getID());
    auto const& solid = medium->phase("Solid");
    auto const& gas = medium->phase("Gas");
    MPL::VariableArray vars;

    auto const T_ref =
        medium->property(MPL::PropertyType::reference_temperature)
            .template value<double>(vars, x_position, t, dt);
    vars[static_cast<int>(MPL::Variable::temperature)] = T_ref;

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
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                _element, N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        auto& eps = _ip_data[ip].eps;
        auto const& sigma_eff = _ip_data[ip].sigma_eff;

        vars[static_cast<int>(MPL::Variable::phase_pressure)] = N_p.dot(p);

        auto const alpha = solid.property(MPL::PropertyType::biot_coefficient)
                               .template value<double>(vars, x_position, t, dt);
        auto const rho_sr =
            solid.property(MPL::PropertyType::density)
                .template value<double>(vars, x_position, t, dt);
        auto const porosity =
            solid.property(MPL::PropertyType::porosity)
                .template value<double>(vars, x_position, t, dt);

        auto const rho_fr =
            gas.property(MPL::PropertyType::density)
                .template value<double>(vars, x_position, t, dt);

        auto const& b = _process_data.specific_body_force;
        auto const& identity2 = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value>::identity2;

        eps.noalias() = B * u;
        vars[static_cast<int>(MaterialPropertyLib::Variable::strain)]
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps);

        auto C = _ip_data[ip].updateConstitutiveRelation(vars, t, x_position,
                                                         dt, u, T_ref);

        local_Jac.noalias() += B.transpose() * C * B * w;

        double p_at_xi = 0.;
        NumLib::shapeFunctionInterpolate(p, N_p, p_at_xi);

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
        const double t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_xdot, const double /*dxdot_dx*/,
        const double /*dx_dx*/, int const process_id,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    // For the equations with pressure
    if (process_id == _process_data.hydraulic_process_id)
    {
        assembleWithJacobianForPressureEquations(t, dt, local_x, local_xdot,
                                                 local_b_data, local_Jac_data);
        return;
    }

    // For the equations with deformation
    assembleWithJacobianForDeformationEquations(t, dt, local_x, local_b_data,
                                                local_Jac_data);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void HydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                  ShapeFunctionPressure, IntegrationMethod,
                                  DisplacementDim>::
    setInitialConditionsConcrete(std::vector<double> const& /*local_x*/,
                                 double const /*t*/,
                                 bool const /*use_monolithic_scheme*/,
                                 int const /*process_id*/)
{
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void HydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                  ShapeFunctionPressure, IntegrationMethod,
                                  DisplacementDim>::
    postNonLinearSolverConcrete(std::vector<double> const& local_x,
                                std::vector<double> const& /*local_xdot*/,
                                double const t, double const dt,
                                bool const use_monolithic_scheme,
                                int const /*process_id*/)
{
    const int displacement_offset =
        use_monolithic_scheme ? displacement_index : 0;

    auto u =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x.data() + displacement_offset,
                                      displacement_size);
    MPL::VariableArray vars;
    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    auto const& medium = _process_data.media_map->getMedium(_element.getID());

    auto const T_ref =
        medium->property(MPL::PropertyType::reference_temperature)
            .template value<double>(MPL::VariableArray(), x_position, t, dt);

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& N_u = _ip_data[ip].N_u;
        auto const& dNdx_u = _ip_data[ip].dNdx_u;

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                _element, N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        auto& eps = _ip_data[ip].eps;
        eps.noalias() = B * u;
        vars[static_cast<int>(MaterialPropertyLib::Variable::strain)]
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps);

        _ip_data[ip].updateConstitutiveRelation(vars, t, x_position, dt, u,
                                                T_ref);
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::size_t HydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::setIPDataInitialConditions(std::string const& name,
                                                 double const* values,
                                                 int const integration_order)
{
    if (integration_order !=
        static_cast<int>(_integration_method.getIntegrationOrder()))
    {
        OGS_FATAL(
            "Setting integration point initial conditions; The integration "
            "order of the local assembler for element {:d} is different from "
            "the integration order in the initial condition.",
            _element.getID());
    }

    if (name == "sigma_ip")
    {
        return setSigma(values);
    }
    if (name == "epsilon_ip")
    {
        return setEpsilon(values);
    }

    return 0;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::size_t HydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::setSigma(double const* values)
{
    return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
        values, _ip_data, &IpData::sigma_eff);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> HydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::getSigma() const
{
    return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
        _ip_data, &IpData::sigma_eff);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::size_t HydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::setEpsilon(double const* values)
{
    return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
        values, _ip_data, &IpData::eps);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> HydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::getEpsilon() const
{
    auto const kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    std::vector<double> ip_epsilon_values;
    auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, Eigen::Dynamic, kelvin_vector_size, Eigen::RowMajor>>(
        ip_epsilon_values, n_integration_points, kelvin_vector_size);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        auto const& eps = _ip_data[ip].eps;
        cache_mat.row(ip) =
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(eps);
    }

    return ip_epsilon_values;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void HydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                  ShapeFunctionPressure, IntegrationMethod,
                                  DisplacementDim>::
    computeSecondaryVariableConcrete(double const t, double const dt,
                                     Eigen::VectorXd const& local_x,
                                     Eigen::VectorXd const& /*local_x_dot*/)
{
    auto const p = local_x.template segment<pressure_size>(pressure_index);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(_element, _is_axially_symmetric, p,
                         *_process_data.pressure_interpolated);

    int const elem_id = _element.getID();
    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(elem_id);
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    auto const& medium = _process_data.media_map->getMedium(elem_id);
    auto const& solid = medium->phase("Solid");
    MPL::VariableArray vars;

    SymmetricTensor k_sum = SymmetricTensor::Zero(KelvinVectorSize);
    auto sigma_eff_sum = MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
        Eigen::Matrix<double, 3, 3>::Zero());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto const& eps = _ip_data[ip].eps;
        auto const& sigma_eff = _ip_data[ip].sigma_eff;
        sigma_eff_sum += sigma_eff;

        auto const alpha_b =
            solid.property(MPL::PropertyType::biot_coefficient)
                .template value<double>(vars, x_position, t, dt);
        auto const& identity2 = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value>::identity2;

        double const p_int_pt = _ip_data[ip].N_p.dot(p);
        vars[static_cast<int>(MPL::Variable::phase_pressure)] = p_int_pt;

        vars[static_cast<int>(MPL::Variable::strain)].emplace<SymmetricTensor>(
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(eps));
        vars[static_cast<int>(MPL::Variable::total_stress)]
            .emplace<SymmetricTensor>(
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                    (sigma_eff - alpha_b * identity2 * p_int_pt).eval()));
        k_sum += MPL::getSymmetricTensor<DisplacementDim>(
                    medium->property(MPL::PropertyType::permeability)
                        .value(vars, x_position, t, dt));
    }

    Eigen::Map<Eigen::VectorXd>(
        &(*_process_data.permeability)[elem_id * KelvinVectorSize],
        KelvinVectorSize) = k_sum / n_integration_points;

    Eigen::Matrix<double, 3, 3, 0, 3, 3> const sigma_avg =
        MathLib::KelvinVector::kelvinVectorToTensor(sigma_eff_sum) /
        n_integration_points;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> e_s(sigma_avg);

    Eigen::Map<Eigen::Vector3d>(
        &(*_process_data.principal_stress_values)[elem_id * 3], 3) =
        e_s.eigenvalues();

    auto eigen_vectors = e_s.eigenvectors();

    for (auto i = 0; i < 3; i++)
    {
        Eigen::Map<Eigen::Vector3d>(
            &(*_process_data.principal_stress_vector[i])[elem_id * 3], 3) =
            eigen_vectors.col(i);
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
unsigned HydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::getNumberOfIntegrationPoints() const
{
    return _integration_method.getNumberOfPoints();
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
typename MaterialLib::Solids::MechanicsBase<
    DisplacementDim>::MaterialStateVariables const&
HydroMechanicsLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                             IntegrationMethod, DisplacementDim>::
    getMaterialStateVariablesAt(unsigned integration_point) const
{
    return *_ip_data[integration_point].material_state_variables;
}

}  // namespace HydroMechanics
}  // namespace ProcessLib
