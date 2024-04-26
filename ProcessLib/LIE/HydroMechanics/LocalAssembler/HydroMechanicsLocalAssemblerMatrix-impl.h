/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "HydroMechanicsLocalAssemblerMatrix.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/KelvinVector.h"
#include "MeshLib/ElementStatus.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/Interpolation.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int GlobalDim>
HydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                   ShapeFunctionPressure, GlobalDim>::
    HydroMechanicsLocalAssemblerMatrix(
        MeshLib::Element const& e,
        std::size_t const n_variables,
        std::size_t const /*local_matrix_size*/,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        HydroMechanicsProcessData<GlobalDim>& process_data)
    : HydroMechanicsLocalAssemblerInterface(
          e, is_axially_symmetric,
          (n_variables - 1) * ShapeFunctionDisplacement::NPOINTS * GlobalDim +
              ShapeFunctionPressure::NPOINTS,
          dofIndex_to_localIndex),
      _process_data(process_data)
{
    unsigned const n_integration_points =
        integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);
    _secondary_data.N.resize(n_integration_points);

    auto const shape_matrices_u =
        NumLib::initShapeMatrices<ShapeFunctionDisplacement,
                                  ShapeMatricesTypeDisplacement, GlobalDim>(
            e, is_axially_symmetric, integration_method);

    auto const shape_matrices_p =
        NumLib::initShapeMatrices<ShapeFunctionPressure,
                                  ShapeMatricesTypePressure, GlobalDim>(
            e, is_axially_symmetric, integration_method);

    auto& solid_material = MaterialLib::Solids::selectSolidConstitutiveRelation(
        _process_data.solid_materials, _process_data.material_ids, e.getID());

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(e.getID());
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        _ip_data.emplace_back(solid_material);
        auto& ip_data = _ip_data[ip];
        auto const& sm_u = shape_matrices_u[ip];
        auto const& sm_p = shape_matrices_p[ip];

        ip_data.integration_weight =
            sm_u.detJ * sm_u.integralMeasure *
            integration_method.getWeightedPoint(ip).getWeight();
        ip_data.darcy_velocity.setZero();

        ip_data.N_u = sm_u.N;
        ip_data.dNdx_u = sm_u.dNdx;
        ip_data.H_u.setZero(GlobalDim, displacement_size);
        for (int i = 0; i < GlobalDim; ++i)
        {
            ip_data.H_u
                .template block<1, displacement_size / GlobalDim>(
                    i, i * displacement_size / GlobalDim)
                .noalias() = ip_data.N_u;
        }

        ip_data.N_p = sm_p.N;
        ip_data.dNdx_p = sm_p.dNdx;

        _secondary_data.N[ip] = sm_u.N;

        ip_data.sigma_eff.setZero(kelvin_vector_size);
        ip_data.eps.setZero(kelvin_vector_size);

        ip_data.sigma_eff_prev.resize(kelvin_vector_size);
        ip_data.eps_prev.resize(kelvin_vector_size);
        ip_data.C.resize(kelvin_vector_size, kelvin_vector_size);

        auto const initial_effective_stress =
            _process_data.initial_effective_stress(0, x_position);
        for (unsigned i = 0; i < kelvin_vector_size; i++)
        {
            ip_data.sigma_eff[i] = initial_effective_stress[i];
            ip_data.sigma_eff_prev[i] = initial_effective_stress[i];
        }
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int GlobalDim>
void HydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                        ShapeFunctionPressure, GlobalDim>::
    assembleWithJacobianConcrete(double const t, double const dt,
                                 Eigen::VectorXd const& local_x,
                                 Eigen::VectorXd const& local_x_prev,
                                 Eigen::VectorXd& local_rhs,
                                 Eigen::MatrixXd& local_Jac)
{
    auto p = const_cast<Eigen::VectorXd&>(local_x).segment(pressure_index,
                                                           pressure_size);
    auto p_prev = const_cast<Eigen::VectorXd&>(local_x_prev)
                      .segment(pressure_index, pressure_size);

    if (_process_data.deactivate_matrix_in_flow)
    {
        setPressureOfInactiveNodes(t, p);
    }

    auto u = local_x.segment(displacement_index, displacement_size);
    auto u_prev = local_x_prev.segment(displacement_index, displacement_size);

    auto rhs_p = local_rhs.template segment<pressure_size>(pressure_index);
    auto rhs_u =
        local_rhs.template segment<displacement_size>(displacement_index);

    auto J_pp = local_Jac.template block<pressure_size, pressure_size>(
        pressure_index, pressure_index);
    auto J_pu = local_Jac.template block<pressure_size, displacement_size>(
        pressure_index, displacement_index);
    auto J_uu = local_Jac.template block<displacement_size, displacement_size>(
        displacement_index, displacement_index);
    auto J_up = local_Jac.template block<displacement_size, pressure_size>(
        displacement_index, pressure_index);

    assembleBlockMatricesWithJacobian(t, dt, p, p_prev, u, u_prev, rhs_p, rhs_u,
                                      J_pp, J_pu, J_uu, J_up);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int GlobalDim>
void HydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                        ShapeFunctionPressure, GlobalDim>::
    assembleBlockMatricesWithJacobian(
        double const t, double const dt,
        Eigen::Ref<const Eigen::VectorXd> const& p,
        Eigen::Ref<const Eigen::VectorXd> const& p_prev,
        Eigen::Ref<const Eigen::VectorXd> const& u,
        Eigen::Ref<const Eigen::VectorXd> const& u_prev,
        Eigen::Ref<Eigen::VectorXd> rhs_p, Eigen::Ref<Eigen::VectorXd> rhs_u,
        Eigen::Ref<Eigen::MatrixXd> J_pp, Eigen::Ref<Eigen::MatrixXd> J_pu,
        Eigen::Ref<Eigen::MatrixXd> J_uu, Eigen::Ref<Eigen::MatrixXd> J_up)
{
    assert(this->_element.getDimension() == GlobalDim);

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

    auto const& gravity_vec = _process_data.specific_body_force;

    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;
    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points = _ip_data.size();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto& ip_data = _ip_data[ip];
        auto const& ip_w = ip_data.integration_weight;
        auto const& N_u = ip_data.N_u;
        auto const& dNdx_u = ip_data.dNdx_u;
        auto const& N_p = ip_data.N_p;
        auto const& dNdx_p = ip_data.dNdx_p;
        auto const& H_u = ip_data.H_u;

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                _element, N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<GlobalDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        auto const& eps_prev = ip_data.eps_prev;
        auto const& sigma_eff_prev = ip_data.sigma_eff_prev;
        auto& sigma_eff = ip_data.sigma_eff;

        auto& eps = ip_data.eps;
        auto& state = ip_data.material_state_variables;

        auto const alpha = _process_data.biot_coefficient(t, x_position)[0];
        auto const rho_sr = _process_data.solid_density(t, x_position)[0];
        auto const rho_fr = _process_data.fluid_density(t, x_position)[0];
        auto const porosity = _process_data.porosity(t, x_position)[0];

        double const rho = rho_sr * (1 - porosity) + porosity * rho_fr;
        auto const& identity2 =
            MathLib::KelvinVector::Invariants<kelvin_vector_size>::identity2;

        eps.noalias() = B * u;

        variables.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<GlobalDim>>(eps);

        variables_prev.stress
            .emplace<MathLib::KelvinVector::KelvinVectorType<GlobalDim>>(
                sigma_eff_prev);
        variables_prev.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<GlobalDim>>(
                eps_prev);
        variables_prev.temperature = _process_data.reference_temperature;

        auto&& solution = _ip_data[ip].solid_material.integrateStress(
            variables_prev, variables, t, x_position, dt, *state);

        if (!solution)
        {
            OGS_FATAL("Computation of local constitutive relation failed.");
        }

        MathLib::KelvinVector::KelvinMatrixType<GlobalDim> C;
        std::tie(sigma_eff, state, C) = std::move(*solution);

        J_uu.noalias() += B.transpose() * C * B * ip_w;

        rhs_u.noalias() -= B.transpose() * sigma_eff * ip_w;
        rhs_u.noalias() -= -H_u.transpose() * rho * gravity_vec * ip_w;

        //
        // pressure equation, pressure part and displacement equation, pressure
        // part
        //
        if (!_process_data.deactivate_matrix_in_flow)  // Only for hydraulically
                                                       // active matrix
        {
            Kup.noalias() += B.transpose() * alpha * identity2 * N_p * ip_w;

            double const k_over_mu =
                _process_data.intrinsic_permeability(t, x_position)[0] /
                _process_data.fluid_viscosity(t, x_position)[0];
            double const S = _process_data.specific_storage(t, x_position)[0];

            auto q = ip_data.darcy_velocity.head(GlobalDim);
            q.noalias() = -k_over_mu * (dNdx_p * p + rho_fr * gravity_vec);

            laplace_p.noalias() +=
                dNdx_p.transpose() * k_over_mu * dNdx_p * ip_w;
            storage_p.noalias() += N_p.transpose() * S * N_p * ip_w;

            rhs_p.noalias() +=
                dNdx_p.transpose() * rho_fr * k_over_mu * gravity_vec * ip_w;
        }
    }

    // displacement equation, pressure part
    J_up.noalias() -= Kup;

    // pressure equation, pressure part.
    J_pp.noalias() += laplace_p + storage_p / dt;

    // pressure equation, displacement part.
    J_pu.noalias() += Kup.transpose() / dt;

    // pressure equation
    rhs_p.noalias() -= laplace_p * p + storage_p * (p - p_prev) / dt +
                       Kup.transpose() * (u - u_prev) / dt;

    // displacement equation
    rhs_u.noalias() -= -Kup * p;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int GlobalDim>
void HydroMechanicsLocalAssemblerMatrix<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    GlobalDim>::postTimestepConcreteWithVector(double const t, double const dt,
                                               Eigen::VectorXd const& local_x)
{
    auto p = const_cast<Eigen::VectorXd&>(local_x).segment(pressure_index,
                                                           pressure_size);
    if (_process_data.deactivate_matrix_in_flow)
    {
        setPressureOfInactiveNodes(t, p);
    }
    auto u = local_x.segment(displacement_index, displacement_size);

    postTimestepConcreteWithBlockVectors(t, dt, p, u);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int GlobalDim>
void HydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                        ShapeFunctionPressure, GlobalDim>::
    postTimestepConcreteWithBlockVectors(
        double const t, double const dt,
        Eigen::Ref<const Eigen::VectorXd> const& p,
        Eigen::Ref<const Eigen::VectorXd> const& u)
{
    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;
    ParameterLib::SpatialPosition x_position;

    auto const e_id = _element.getID();
    x_position.setElementID(e_id);

    using KV = MathLib::KelvinVector::KelvinVectorType<GlobalDim>;
    KV sigma_avg = KV::Zero();
    GlobalDimVector velocity_avg;
    velocity_avg.setZero();

    unsigned const n_integration_points = _ip_data.size();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto& ip_data = _ip_data[ip];

        auto const& eps_prev = ip_data.eps_prev;
        auto const& sigma_eff_prev = ip_data.sigma_eff_prev;

        auto& eps = ip_data.eps;
        auto& sigma_eff = ip_data.sigma_eff;
        auto& state = ip_data.material_state_variables;

        auto const& N_u = ip_data.N_u;
        auto const& dNdx_u = ip_data.dNdx_u;

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                _element, N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<GlobalDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        eps.noalias() = B * u;

        variables.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<GlobalDim>>(eps);

        variables_prev.stress
            .emplace<MathLib::KelvinVector::KelvinVectorType<GlobalDim>>(
                sigma_eff_prev);
        variables_prev.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<GlobalDim>>(
                eps_prev);
        variables_prev.temperature = _process_data.reference_temperature;

        auto&& solution = _ip_data[ip].solid_material.integrateStress(
            variables_prev, variables, t, x_position, dt, *state);

        if (!solution)
        {
            OGS_FATAL("Computation of local constitutive relation failed.");
        }

        MathLib::KelvinVector::KelvinMatrixType<GlobalDim> C;
        std::tie(sigma_eff, state, C) = std::move(*solution);

        sigma_avg += ip_data.sigma_eff;

        if (!_process_data.deactivate_matrix_in_flow)  // Only for hydraulically
                                                       // active matrix
        {
            double const k_over_mu =
                _process_data.intrinsic_permeability(t, x_position)[0] /
                _process_data.fluid_viscosity(t, x_position)[0];
            auto const rho_fr = _process_data.fluid_density(t, x_position)[0];
            auto const& gravity_vec = _process_data.specific_body_force;
            auto const& dNdx_p = ip_data.dNdx_p;

            ip_data.darcy_velocity.head(GlobalDim).noalias() =
                -k_over_mu * (dNdx_p * p + rho_fr * gravity_vec);
            velocity_avg += ip_data.darcy_velocity.head(GlobalDim);
        }
    }

    sigma_avg /= n_integration_points;
    velocity_avg /= n_integration_points;

    Eigen::Map<KV>(
        &(*_process_data.element_stresses)[e_id * KV::RowsAtCompileTime]) =
        MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_avg);

    Eigen::Map<GlobalDimVector>(
        &(*_process_data.element_velocities)[e_id * GlobalDim]) = velocity_avg;

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        GlobalDim>(_element, _is_axially_symmetric, p,
                   *_process_data.mesh_prop_nodal_p);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int GlobalDim>
void HydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                        ShapeFunctionPressure, GlobalDim>::
    setPressureOfInactiveNodes(double const t, Eigen::Ref<Eigen::VectorXd> p)
{
    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());
    for (unsigned i = 0; i < pressure_size; i++)
    {
        // only inactive nodes
        if (_process_data.p_element_status->isActiveNode(_element.getNode(i)))
        {
            continue;
        }
        x_position.setNodeID(getNodeIndex(_element, i));
        auto const p0 = (*_process_data.p0)(t, x_position)[0];
        p[i] = p0;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int GlobalDim>
std::vector<double> const& HydroMechanicsLocalAssemblerMatrix<
    ShapeFunctionDisplacement, ShapeFunctionPressure, GlobalDim>::
    getIntPtSigma(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointKelvinVectorData<GlobalDim>(
        _ip_data, &IntegrationPointDataType::sigma_eff, cache);
}
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int GlobalDim>
std::vector<double> const& HydroMechanicsLocalAssemblerMatrix<
    ShapeFunctionDisplacement, ShapeFunctionPressure, GlobalDim>::
    getIntPtEpsilon(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointKelvinVectorData<GlobalDim>(
        _ip_data, &IntegrationPointDataType::eps, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int GlobalDim>
std::vector<double> const& HydroMechanicsLocalAssemblerMatrix<
    ShapeFunctionDisplacement, ShapeFunctionPressure, GlobalDim>::
    getIntPtDarcyVelocity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    unsigned const n_integration_points = _ip_data.size();

    cache.clear();
    auto cache_matrix = MathLib::createZeroedMatrix<
        Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, GlobalDim, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        cache_matrix.col(ip).noalias() = _ip_data[ip].darcy_velocity;
    }

    return cache;
}
}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
