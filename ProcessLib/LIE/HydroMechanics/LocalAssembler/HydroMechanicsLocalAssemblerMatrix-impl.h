/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_LIE_HYDROMECHANICS_HYDROMECHANICSLOCALASSEMBLER_MATRIX_IMPL_H_
#define PROCESSLIB_LIE_HYDROMECHANICS_HYDROMECHANICSLOCALASSEMBLER_MATRIX_IMPL_H_

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
          typename IntegrationMethod, unsigned GlobalDim>
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
          e,
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
        initShapeMatrices<ShapeFunctionDisplacement, ShapeMatricesTypeDisplacement,
                          IntegrationMethod, GlobalDim>(
                e, is_axially_symmetric, integration_method);

    auto const shape_matrices_p =
        initShapeMatrices<ShapeFunctionPressure, ShapeMatricesTypePressure,
                          IntegrationMethod, GlobalDim>(
                e, is_axially_symmetric, integration_method);

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
        ip_data.b_matrices.resize(kelvin_vector_size, displacement_size);

        ip_data.N_u = sm_u.N;
        ip_data.H_u.setZero(GlobalDim, displacement_size);
        for (unsigned i = 0; i < GlobalDim; ++i)
            ip_data.H_u.template block<1, displacement_size / GlobalDim>(
                   i, i * displacement_size / GlobalDim)
                .noalias() = ip_data.N_u;
        auto const x_coord =
            interpolateXCoordinate<ShapeFunctionDisplacement,
                                   ShapeMatricesTypeDisplacement>(e, sm_u.N);
        LinearBMatrix::computeBMatrix<GlobalDim,
                                      ShapeFunctionDisplacement::NPOINTS>(
            sm_u.dNdx, ip_data.b_matrices, is_axially_symmetric, sm_u.N,
            x_coord);

        ip_data.N_p = sm_p.N;
        ip_data.dNdx_p = sm_p.dNdx;

        ip_data.sigma_eff.resize(kelvin_vector_size);
        ip_data.sigma_eff_prev.resize(kelvin_vector_size);
        ip_data.eps.resize(kelvin_vector_size);
        ip_data.eps_prev.resize(kelvin_vector_size);
        ip_data.C.resize(kelvin_vector_size, kelvin_vector_size);

        auto const initial_effective_stress = _process_data.initial_effective_stress(0, x_position);
        for (unsigned i=0; i<kelvin_vector_size; i++)
        {
            ip_data.sigma_eff[i] = initial_effective_stress[i];
            ip_data.sigma_eff_prev[i] = initial_effective_stress[i];
        }
    }

    if (_process_data.use_initial_stress_as_reference)
    {
        _initial_pressure.resize(pressure_size);
        for (unsigned i=0; i<pressure_size; i++)
        {
            x_position.setNodeID(e.getNodeIndex(i));
            _initial_pressure[i] = _process_data.initial_pressure(0, x_position)[0];
        }
    }
}


template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, unsigned GlobalDim>
void
HydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                   ShapeFunctionPressure, IntegrationMethod,
                                   GlobalDim>::
assembleWithJacobianConcrete(
    double const t,
    Eigen::VectorXd const& local_x,
    Eigen::VectorXd const& local_x_dot,
    Eigen::VectorXd& local_rhs,
    Eigen::MatrixXd& local_Jac)
{
    auto p = const_cast<Eigen::VectorXd&>(local_x).segment(pressure_index, pressure_size);
    auto p_dot = const_cast<Eigen::VectorXd&>(local_x_dot).segment(pressure_index, pressure_size);

    if (_process_data.pv_p && !_process_data.pv_p->getElementStatus().isActive(_element.getID()))
    {
        SpatialPosition x_position;
        x_position.setElementID(_element.getID());
        for (unsigned i=0; i<pressure_size; i++)
        {
            // only inactive nodes
            if (_process_data.pv_p->getElementStatus().isActiveNode(_element.getNode(i)))
                continue;
            x_position.setNodeID(_element.getNodeIndex(i));
            auto const p0 = _process_data.pv_p->getInitialCondition()(t, x_position)[0];
            p[i] = p0;
            p_dot[i] = 0;
        }
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
          typename IntegrationMethod, unsigned GlobalDim>
void
HydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                   ShapeFunctionPressure, IntegrationMethod,
                                   GlobalDim>::
assembleBlockMatricesWithJacobian(
    double const t,
    Eigen::Ref<const Eigen::VectorXd> const& p,
    Eigen::Ref<const Eigen::VectorXd> const& p_dot,
    Eigen::Ref<const Eigen::VectorXd> const& u,
    Eigen::Ref<const Eigen::VectorXd> const& u_dot,
    Eigen::Ref<Eigen::VectorXd> rhs_p,
    Eigen::Ref<Eigen::VectorXd> rhs_u,
    Eigen::Ref<Eigen::MatrixXd> J_pp,
    Eigen::Ref<Eigen::MatrixXd> J_pu,
    Eigen::Ref<Eigen::MatrixXd> J_uu,
    Eigen::Ref<Eigen::MatrixXd> J_up)
{
    assert (this->_element.getDimension() == GlobalDim);

    typename ShapeMatricesTypePressure::NodalMatrixType laplace_p;
    laplace_p.setZero(pressure_size, pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage_p;
    storage_p.setZero(pressure_size, pressure_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        displacement_size, pressure_size>
        Kup;
    Kup.setZero(displacement_size, pressure_size);

    double const& dt = _process_data.dt;

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points = _ip_data.size();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto& ip_data = _ip_data[ip];
        auto const& ip_w = ip_data.integration_weight;
        auto const& N_p = ip_data.N_p;
        auto const& dNdx_p = ip_data.dNdx_p;
        auto const& H_u = ip_data.H_u;

        auto const& B = ip_data.b_matrices;
        auto const& eps_prev = ip_data.eps_prev;
        auto const& sigma_eff_prev = ip_data.sigma_eff_prev;

        auto& eps = ip_data.eps;
        auto& sigma_eff = ip_data.sigma_eff;
        auto& C = ip_data.C;
        auto& material_state_variables = *ip_data.material_state_variables;

        auto q = ip_data.darcy_velocity.head(GlobalDim);

        double const S =
            _process_data.specific_storage(t, x_position)[0];
        double const k_over_mu =
            _process_data.intrinsic_permeability(t, x_position)[0] /
            _process_data.fluid_viscosity(t, x_position)[0];
        auto const alpha = _process_data.biot_coefficient(t, x_position)[0];
        auto const rho_sr = _process_data.solid_density(t, x_position)[0];
        auto const rho_fr = _process_data.fluid_density(t, x_position)[0];
        auto const porosity = _process_data.porosity(t, x_position)[0];
        auto const& gravity_vec = _process_data.specific_body_force;

        double const rho = rho_sr * (1 - porosity) + porosity * rho_fr;
        auto const& identity2 =
            MaterialLib::SolidModels::Invariants<kelvin_vector_size>::identity2;

        eps.noalias() = B * u;

        if (!_ip_data[ip].solid_material.computeConstitutiveRelation(
                t, x_position, _process_data.dt, eps_prev, eps, sigma_eff_prev,
                sigma_eff, C, material_state_variables))
            OGS_FATAL("Computation of local constitutive relation failed.");

        q.noalias() = - k_over_mu * (dNdx_p * p + rho_fr * gravity_vec);

        J_uu.noalias() += B.transpose() * C * B * ip_w;

        if (_process_data.use_initial_stress_as_reference)
        {
            auto const vec_sigma_eff_ref = _process_data.initial_effective_stress(t, x_position);
            auto const sigma_eff_ref =
                Eigen::Map<typename BMatricesType::KelvinVectorType const>(vec_sigma_eff_ref.data(), kelvin_vector_size);
            rhs_u.noalias() -= B.transpose() * (sigma_eff - sigma_eff_ref) * ip_w;
        }
        else
        {
            rhs_u.noalias() -= B.transpose() * sigma_eff * ip_w;
        }
        rhs_u.noalias() -= - H_u.transpose() * rho * gravity_vec * ip_w;

        //
        // displacement equation, pressure part
        //
        Kup.noalias() += B.transpose() * alpha * identity2 * N_p * ip_w;

        //
        // pressure equation, pressure part.
        //
        laplace_p.noalias() += dNdx_p.transpose() * k_over_mu * dNdx_p * ip_w;

        storage_p.noalias() += N_p.transpose() * S * N_p * ip_w;

        rhs_p.noalias() += dNdx_p.transpose() * rho_fr * k_over_mu * gravity_vec * ip_w;
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
    if (_process_data.use_initial_stress_as_reference)
        rhs_u.noalias() -= - Kup * (p - _initial_pressure);
    else
        rhs_u.noalias() -= - Kup * p;
}


template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, unsigned GlobalDim>
void
HydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                   ShapeFunctionPressure, IntegrationMethod,
                                   GlobalDim>::
computeSecondaryVariableConcreteWithVector(
    double const t,
    Eigen::VectorXd const& local_x)
{
    auto p = const_cast<Eigen::VectorXd&>(local_x).segment(pressure_index, pressure_size);

    if (_process_data.pv_p && !_process_data.pv_p->getElementStatus().isActive(_element.getID()))
    {
        SpatialPosition x_position;
        x_position.setElementID(_element.getID());
        for (unsigned i=0; i<pressure_size; i++)
        {
            // only inactive nodes
            if (_process_data.pv_p->getElementStatus().isActiveNode(_element.getNode(i)))
                continue;
            x_position.setNodeID(_element.getNodeIndex(i));
            auto const p0 = _process_data.pv_p->getInitialCondition()(t, x_position)[0];
            p[i] = p0;
        }
    }

    auto u = local_x.segment(displacement_index, displacement_size);

    computeSecondaryVariableConcreteWithBlockVectors(t, p, u);
}



template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, unsigned GlobalDim>
void
HydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                   ShapeFunctionPressure, IntegrationMethod,
                                   GlobalDim>::
computeSecondaryVariableConcreteWithBlockVectors(
    double const t,
    Eigen::Ref<const Eigen::VectorXd> const& /*p*/,
    Eigen::Ref<const Eigen::VectorXd> const& u)
{
    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points = _ip_data.size();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto& ip_data = _ip_data[ip];

        auto const& B = ip_data.b_matrices;
        auto const& eps_prev = ip_data.eps_prev;
        auto const& sigma_eff_prev = ip_data.sigma_eff_prev;

        auto& eps = ip_data.eps;
        auto& sigma_eff = ip_data.sigma_eff;
        auto& C = ip_data.C;
        auto& material_state_variables = *ip_data.material_state_variables;

        //auto q = ip_data.darcy_velocity.head(GlobalDim);

        eps.noalias() = B * u;

        if (!_ip_data[ip].solid_material.computeConstitutiveRelation(
                t, x_position, _process_data.dt, eps_prev, eps, sigma_eff_prev,
                sigma_eff, C, material_state_variables))
            OGS_FATAL("Computation of local constitutive relation failed.");

        //q.noalias() = - k_over_mu * (dNdx_p * p + rho_fr * gravity_vec);
    }



    Eigen::Vector3d ele_stress;
    ele_stress.setZero();
    Eigen::Vector3d ele_strain;
    ele_strain.setZero();
//    Eigen::Vector3d ele_velocity;
//    ele_velocity.setZero();

    for (auto const& ip_data : _ip_data)
    {
        ele_stress[0] += ip_data.sigma_eff[0];
        ele_stress[1] += ip_data.sigma_eff[1];
        ele_stress[2] += ip_data.sigma_eff[3];

        ele_strain[0] += ip_data.eps[0];
        ele_strain[1] += ip_data.eps[1];
        ele_strain[2] += ip_data.eps[3];

//        ele_velocity += ip_data.darcy_velocity;
    }

    ele_stress /= _ip_data.size();
    ele_strain /= _ip_data.size();
//    ele_velocity /= _ip_data.size();

    (*_process_data.mesh_prop_stress_xx)[_element.getID()] = ele_stress[0];
    (*_process_data.mesh_prop_stress_yy)[_element.getID()] = ele_stress[1];
    (*_process_data.mesh_prop_stress_xy)[_element.getID()] = ele_stress[2];
    (*_process_data.mesh_prop_strain_xx)[_element.getID()] = ele_strain[0];
    (*_process_data.mesh_prop_strain_yy)[_element.getID()] = ele_strain[1];
    (*_process_data.mesh_prop_strain_xy)[_element.getID()] = ele_strain[2];
//    for (unsigned i=0; i<3; i++)
//        (*_process_data.mesh_prop_velocity)[_element.getID()*3 + i] = ele_velocity[i];
}


template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, unsigned GlobalDim>
void
HydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                   ShapeFunctionPressure, IntegrationMethod,
                                   GlobalDim>::
postTimestepConcrete(std::vector<double> const& /*local_x*/)
{
    Eigen::Vector3d ele_stress;
    ele_stress.setZero();
    Eigen::Vector3d ele_strain;
    ele_strain.setZero();
    Eigen::Vector3d ele_velocity;
    ele_velocity.setZero();

    for (auto const& ip_data : _ip_data)
    {
        ele_stress[0] += ip_data.sigma_eff[0];
        ele_stress[1] += ip_data.sigma_eff[1];
        // sigma_xy
        ele_stress[2] += ip_data.sigma_eff[3]; // sigma_eff[2] stores sigma_z

        ele_strain[0] += ip_data.eps[0];
        ele_strain[1] += ip_data.eps[1];
        // strain_xy
        ele_strain[2] += ip_data.eps[3];

        ele_velocity += ip_data.darcy_velocity;
    }

    ele_stress /= _ip_data.size();
    ele_strain /= _ip_data.size();
    ele_velocity /= _ip_data.size();

    (*_process_data.mesh_prop_stress_xx)[_element.getID()] = ele_stress[0];
    (*_process_data.mesh_prop_stress_yy)[_element.getID()] = ele_stress[1];
    (*_process_data.mesh_prop_stress_xy)[_element.getID()] = ele_stress[2];
    (*_process_data.mesh_prop_strain_xx)[_element.getID()] = ele_strain[0];
    (*_process_data.mesh_prop_strain_yy)[_element.getID()] = ele_strain[1];
    (*_process_data.mesh_prop_strain_xy)[_element.getID()] = ele_strain[2];
    for (unsigned i=0; i<3; i++)
        (*_process_data.mesh_prop_velocity)[_element.getID()*3 + i] = ele_velocity[i];
}

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib

#endif  // PROCESSLIB_LIE_HYDROMECHANICS_HYDROMECHANICSLOCALASSEMBLER_MATRIX_IMPL_H_
