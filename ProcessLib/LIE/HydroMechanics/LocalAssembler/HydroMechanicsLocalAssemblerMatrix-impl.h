/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "HydroMechanicsLocalAssemblerMatrix.h"

#include "MaterialLib/SolidModels/KelvinVector.h"
#include "MeshLib/ElementStatus.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int GlobalDim>
HydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                   ShapeFunctionPressure, IntegrationMethod,
                                   GlobalDim>::
    HydroMechanicsLocalAssemblerMatrix(
        MeshLib::Element const& e,
        std::size_t const n_variables,
        std::size_t const /*local_matrix_size*/,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        HydroMechanicsProcessData<GlobalDim>& process_data)
    : HydroMechanicsLocalAssemblerInterface(
          e, is_axially_symmetric,
          (n_variables - 1) * ShapeFunctionDisplacement::NPOINTS * GlobalDim +
              ShapeFunctionPressure::NPOINTS,
          dofIndex_to_localIndex),
      _process_data(process_data)
{
    IntegrationMethod integration_method(integration_order);
    unsigned const n_integration_points =
        integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);

    auto const shape_matrices_u =
        initShapeMatrices<ShapeFunctionDisplacement,
                          ShapeMatricesTypeDisplacement, IntegrationMethod,
                          GlobalDim>(e, is_axially_symmetric,
                                     integration_method);

    auto const shape_matrices_p =
        initShapeMatrices<ShapeFunctionPressure, ShapeMatricesTypePressure,
                          IntegrationMethod, GlobalDim>(e, is_axially_symmetric,
                                                        integration_method);

    SpatialPosition x_position;
    x_position.setElementID(e.getID());
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        _ip_data.emplace_back(*_process_data.material);
        auto& ip_data = _ip_data[ip];
        auto const& sm_u = shape_matrices_u[ip];
        auto const& sm_p = shape_matrices_p[ip];
        ip_data.integration_weight =
            sm_u.detJ * sm_u.integralMeasure *
            integration_method.getWeightedPoint(ip).getWeight();

        ip_data.N_u = sm_u.N;
        ip_data.dNdx_u = sm_u.dNdx;
        ip_data.H_u.setZero(GlobalDim, displacement_size);
        for (int i = 0; i < GlobalDim; ++i)
            ip_data.H_u
                .template block<1, displacement_size / GlobalDim>(
                    i, i * displacement_size / GlobalDim)
                .noalias() = ip_data.N_u;

        ip_data.N_p = sm_p.N;
        ip_data.dNdx_p = sm_p.dNdx;

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
          typename IntegrationMethod, int GlobalDim>
void HydroMechanicsLocalAssemblerMatrix<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    GlobalDim>::assembleWithJacobianConcrete(double const t,
                                             Eigen::VectorXd const& local_x,
                                             Eigen::VectorXd const& local_x_dot,
                                             Eigen::VectorXd& local_rhs,
                                             Eigen::MatrixXd& local_Jac)
{
    auto p = const_cast<Eigen::VectorXd&>(local_x).segment(pressure_index,
                                                           pressure_size);
    auto p_dot = const_cast<Eigen::VectorXd&>(local_x_dot)
                     .segment(pressure_index, pressure_size);

    if (_process_data.deactivate_matrix_in_flow)
    {
        setPressureOfInactiveNodes(t, p);
        setPressureDotOfInactiveNodes(p_dot);
    }

    auto u = local_x.segment(displacement_index, displacement_size);
    auto u_dot = local_x_dot.segment(displacement_index, displacement_size);

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

    assembleBlockMatricesWithJacobian(t, p, p_dot, u, u_dot, rhs_p, rhs_u, J_pp,
                                      J_pu, J_uu, J_up);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int GlobalDim>
void HydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                        ShapeFunctionPressure,
                                        IntegrationMethod, GlobalDim>::
    assembleBlockMatricesWithJacobian(
        double const t,
        Eigen::Ref<const Eigen::VectorXd> const& p,
        Eigen::Ref<const Eigen::VectorXd> const& p_dot,
        Eigen::Ref<const Eigen::VectorXd> const& u,
        Eigen::Ref<const Eigen::VectorXd> const& u_dot,
        Eigen::Ref<Eigen::VectorXd>
            rhs_p,
        Eigen::Ref<Eigen::VectorXd>
            rhs_u,
        Eigen::Ref<Eigen::MatrixXd>
            J_pp,
        Eigen::Ref<Eigen::MatrixXd>
            J_pu,
        Eigen::Ref<Eigen::MatrixXd>
            J_uu,
        Eigen::Ref<Eigen::MatrixXd>
            J_up)
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

    double const& dt = _process_data.dt;
    auto const& gravity_vec = _process_data.specific_body_force;

    SpatialPosition x_position;
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
            interpolateXCoordinate<ShapeFunctionDisplacement,
                                   ShapeMatricesTypeDisplacement>(_element,
                                                                  N_u);
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

        auto q = ip_data.darcy_velocity.head(GlobalDim);

        double const S = _process_data.specific_storage(t, x_position)[0];
        double const k_over_mu =
            _process_data.intrinsic_permeability(t, x_position)[0] /
            _process_data.fluid_viscosity(t, x_position)[0];
        auto const alpha = _process_data.biot_coefficient(t, x_position)[0];
        auto const rho_sr = _process_data.solid_density(t, x_position)[0];
        auto const rho_fr = _process_data.fluid_density(t, x_position)[0];
        auto const porosity = _process_data.porosity(t, x_position)[0];

        double const rho = rho_sr * (1 - porosity) + porosity * rho_fr;
        auto const& identity2 =
            MaterialLib::SolidModels::Invariants<kelvin_vector_size>::identity2;

        eps.noalias() = B * u;

        auto&& solution = _ip_data[ip].solid_material.integrateStress(
            t, x_position, _process_data.dt, eps_prev, eps, sigma_eff_prev,
            *state);

        if (!solution)
            OGS_FATAL("Computation of local constitutive relation failed.");

        KelvinMatrixType<GlobalDim> C;
        std::tie(sigma_eff, state, C) = std::move(*solution);

        q.noalias() = -k_over_mu * (dNdx_p * p + rho_fr * gravity_vec);

        J_uu.noalias() += B.transpose() * C * B * ip_w;

        rhs_u.noalias() -= B.transpose() * sigma_eff * ip_w;
        rhs_u.noalias() -= -H_u.transpose() * rho * gravity_vec * ip_w;

        //
        // displacement equation, pressure part
        //
        Kup.noalias() += B.transpose() * alpha * identity2 * N_p * ip_w;

        //
        // pressure equation, pressure part.
        //
        laplace_p.noalias() += dNdx_p.transpose() * k_over_mu * dNdx_p * ip_w;

        storage_p.noalias() += N_p.transpose() * S * N_p * ip_w;

        rhs_p.noalias() +=
            dNdx_p.transpose() * rho_fr * k_over_mu * gravity_vec * ip_w;
    }

    // displacement equation, pressure part
    J_up.noalias() -= Kup;

    // pressure equation, pressure part.
    J_pp.noalias() += laplace_p + storage_p / dt;

    // pressure equation, displacement part.
    J_pu.noalias() += Kup.transpose() / dt;

    // pressure equation
    rhs_p.noalias() -=
        laplace_p * p + storage_p * p_dot + Kup.transpose() * u_dot;

    // displacement equation
    rhs_u.noalias() -= -Kup * p;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int GlobalDim>
void HydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                        ShapeFunctionPressure,
                                        IntegrationMethod, GlobalDim>::
    computeSecondaryVariableConcreteWithVector(double const t,
                                               Eigen::VectorXd const& local_x)
{
    auto p = const_cast<Eigen::VectorXd&>(local_x).segment(pressure_index,
                                                           pressure_size);
    if (_process_data.deactivate_matrix_in_flow)
        setPressureOfInactiveNodes(t, p);
    auto u = local_x.segment(displacement_index, displacement_size);

    computeSecondaryVariableConcreteWithBlockVectors(t, p, u);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int GlobalDim>
void HydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                        ShapeFunctionPressure,
                                        IntegrationMethod, GlobalDim>::
    computeSecondaryVariableConcreteWithBlockVectors(
        double const t,
        Eigen::Ref<const Eigen::VectorXd> const& p,
        Eigen::Ref<const Eigen::VectorXd> const& u)
{
    auto const& gravity_vec = _process_data.specific_body_force;

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points = _ip_data.size();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto& ip_data = _ip_data[ip];

        auto const& dNdx_p = ip_data.dNdx_p;
        auto const& eps_prev = ip_data.eps_prev;
        auto const& sigma_eff_prev = ip_data.sigma_eff_prev;

        auto& eps = ip_data.eps;
        auto& sigma_eff = ip_data.sigma_eff;
        auto& state = ip_data.material_state_variables;
        double const k_over_mu =
            _process_data.intrinsic_permeability(t, x_position)[0] /
            _process_data.fluid_viscosity(t, x_position)[0];
        auto const rho_fr = _process_data.fluid_density(t, x_position)[0];

        auto q = ip_data.darcy_velocity.head(GlobalDim);

        auto const& N_u = ip_data.N_u;
        auto const& dNdx_u = ip_data.dNdx_u;

        auto const x_coord =
            interpolateXCoordinate<ShapeFunctionDisplacement,
                                   ShapeMatricesTypeDisplacement>(_element,
                                                                  N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<GlobalDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        eps.noalias() = B * u;

        auto&& solution = _ip_data[ip].solid_material.integrateStress(
            t, x_position, _process_data.dt, eps_prev, eps, sigma_eff_prev,
            *state);

        if (!solution)
            OGS_FATAL("Computation of local constitutive relation failed.");

        KelvinMatrixType<GlobalDim> C;
        std::tie(sigma_eff, state, C) = std::move(*solution);

        q.noalias() = -k_over_mu * (dNdx_p * p + rho_fr * gravity_vec);
    }

    int n = GlobalDim == 2 ? 4 : 6;
    Eigen::VectorXd ele_stress = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd ele_strain = Eigen::VectorXd::Zero(n);
    Eigen::Vector3d ele_velocity = Eigen::Vector3d::Zero();

    for (auto const& ip_data : _ip_data)
    {
        ele_stress += ip_data.sigma_eff;
        ele_strain += ip_data.eps;
        ele_velocity += ip_data.darcy_velocity;
    }

    ele_stress /= static_cast<double>(n_integration_points);
    ele_strain /= static_cast<double>(n_integration_points);
    ele_velocity /= static_cast<double>(n_integration_points);

    auto const element_id = _element.getID();
    (*_process_data.mesh_prop_stress_xx)[_element.getID()] = ele_stress[0];
    (*_process_data.mesh_prop_stress_yy)[_element.getID()] = ele_stress[1];
    (*_process_data.mesh_prop_stress_zz)[_element.getID()] = ele_stress[2];
    (*_process_data.mesh_prop_stress_xy)[_element.getID()] = ele_stress[3];
    if (GlobalDim == 3)
    {
        (*_process_data.mesh_prop_stress_yz)[_element.getID()] = ele_stress[4];
        (*_process_data.mesh_prop_stress_xz)[_element.getID()] = ele_stress[5];
    }

    (*_process_data.mesh_prop_strain_xx)[_element.getID()] = ele_strain[0];
    (*_process_data.mesh_prop_strain_yy)[_element.getID()] = ele_strain[1];
    (*_process_data.mesh_prop_strain_zz)[_element.getID()] = ele_strain[2];
    (*_process_data.mesh_prop_strain_xy)[_element.getID()] = ele_strain[3];
    if (GlobalDim == 3)
    {
        (*_process_data.mesh_prop_strain_yz)[_element.getID()] = ele_strain[4];
        (*_process_data.mesh_prop_strain_xz)[_element.getID()] = ele_strain[5];
    }

    for (unsigned i = 0; i < 3; i++)
        (*_process_data.mesh_prop_velocity)[element_id * 3 + i] =
            ele_velocity[i];
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int GlobalDim>
void HydroMechanicsLocalAssemblerMatrix<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    GlobalDim>::setPressureOfInactiveNodes(double const t,
                                           Eigen::Ref<Eigen::VectorXd>
                                               p)
{
    SpatialPosition x_position;
    x_position.setElementID(_element.getID());
    for (unsigned i = 0; i < pressure_size; i++)
    {
        // only inactive nodes
        if (_process_data.p_element_status->isActiveNode(_element.getNode(i)))
            continue;
        x_position.setNodeID(_element.getNodeIndex(i));
        auto const p0 = (*_process_data.p0)(t, x_position)[0];
        p[i] = p0;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int GlobalDim>
void HydroMechanicsLocalAssemblerMatrix<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    GlobalDim>::setPressureDotOfInactiveNodes(Eigen::Ref<Eigen::VectorXd> p_dot)
{
    for (unsigned i = 0; i < pressure_size; i++)
    {
        // only inactive nodes
        if (_process_data.p_element_status->isActiveNode(_element.getNode(i)))
            continue;
        p_dot[i] = 0;
    }
}

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
