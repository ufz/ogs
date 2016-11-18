/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_LIE_HYDROMECHANICS_HYDROMECHANICSLOCALASSEMBLER_FRACTURE_IMPL_H_
#define PROCESSLIB_LIE_HYDROMECHANICS_HYDROMECHANICSLOCALASSEMBLER_FRACTURE_IMPL_H_

#include "HydroMechanicsLocalAssemblerFracture.h"

#include <iostream>

#include "MaterialLib/FractureModels/FractureIdentity2.h"

#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "ProcessLib/LIE/Common/LevelSetFunction.h"


namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, unsigned GlobalDim>
HydroMechanicsLocalAssemblerFracture<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, IntegrationMethod,
                                     GlobalDim>::
    HydroMechanicsLocalAssemblerFracture(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        HydroMechanicsProcessData<GlobalDim>& process_data)
    : HydroMechanicsLocalAssemblerInterface(
          e,
          ShapeFunctionDisplacement::NPOINTS * GlobalDim + ShapeFunctionPressure::NPOINTS,
          dofIndex_to_localIndex),
      _process_data(process_data)
{
    assert(e.getDimension() == GlobalDim-1);

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

    auto const& frac_prop = *_process_data.fracture_property.get();

    SpatialPosition x_position;
    x_position.setElementID(e.getID());
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        _ip_data.emplace_back(*_process_data.fracture_model);
        auto const& sm_u = shape_matrices_u[ip];
        auto const& sm_p = shape_matrices_p[ip];
        auto& ip_data = _ip_data[ip];
        ip_data.integration_weight =
            sm_u.detJ * sm_u.integralMeasure *
            integration_method.getWeightedPoint(ip).getWeight();

        ip_data.H_u.resize(GlobalDim,
                           ShapeFunctionDisplacement::NPOINTS * GlobalDim);
        computeHMatrix<
            GlobalDim, ShapeFunctionDisplacement::NPOINTS,
            typename ShapeMatricesTypeDisplacement::NodalRowVectorType, HMatrixType>(
            sm_u.N, ip_data.H_u);
        ip_data.N_p = sm_p.N;
        ip_data.dNdx_p = sm_p.dNdx;

        ip_data.w.resize(GlobalDim);
        ip_data.w_prev.resize(GlobalDim);
        ip_data.sigma_eff.resize(GlobalDim);
        ip_data.sigma_eff_prev.resize(GlobalDim);
        ip_data.C.resize(GlobalDim, GlobalDim);
        ip_data.aperture0 = (*frac_prop.aperture0)(0, x_position)[0];
        ip_data.aperture = ip_data.aperture0;

        auto const initial_effective_stress = _process_data.initial_fracture_effective_stress(0, x_position);
        for (unsigned i=0; i<GlobalDim; i++)
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
HydroMechanicsLocalAssemblerFracture<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, IntegrationMethod,
                                     GlobalDim>::
assembleWithJacobianConcrete(
    double const t,
    Eigen::VectorXd const& local_x,
    Eigen::VectorXd const& local_xdot,
    Eigen::VectorXd& local_b,
    Eigen::MatrixXd& local_J)
{
    auto const p = local_x.segment(pressure_index, pressure_size);
    auto const p_dot = local_xdot.segment(pressure_index, pressure_size);
    auto const g = local_x.segment(displacement_index, displacement_size);
    auto const g_dot = local_xdot.segment(displacement_index, displacement_size);

    auto rhs_p = local_b.segment(pressure_index, pressure_size);
    auto rhs_g = local_b.segment(displacement_index, displacement_size);
    auto J_pp = local_J.block(pressure_index, pressure_index, pressure_size,
                              pressure_size);
    auto J_pg = local_J.block(pressure_index, displacement_index, pressure_size,
                              displacement_size);
    auto J_gp = local_J.block(displacement_index, pressure_index,
                              displacement_size, pressure_size);
    auto J_gg = local_J.block(displacement_index, displacement_index,
                              displacement_size, displacement_size);

    assembleBlockMatricesWithJacobian(t, p, p_dot, g, g_dot, rhs_p, rhs_g, J_pp, J_pg, J_gg, J_gp);
}


template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, unsigned GlobalDim>
void HydroMechanicsLocalAssemblerFracture<ShapeFunctionDisplacement,
                                          ShapeFunctionPressure,
                                          IntegrationMethod, GlobalDim>::
assembleBlockMatricesWithJacobian(
    double const t,
    Eigen::Ref<const Eigen::VectorXd> const& nodal_p,
    Eigen::Ref<const Eigen::VectorXd> const& nodal_p_dot,
    Eigen::Ref<const Eigen::VectorXd> const& nodal_g,
    Eigen::Ref<const Eigen::VectorXd> const& nodal_g_dot,
    Eigen::Ref<Eigen::VectorXd> rhs_p,
    Eigen::Ref<Eigen::VectorXd> rhs_g,
    Eigen::Ref<Eigen::MatrixXd> J_pp,
    Eigen::Ref<Eigen::MatrixXd> J_pg,
    Eigen::Ref<Eigen::MatrixXd> J_gg,
    Eigen::Ref<Eigen::MatrixXd> J_gp)
{
    FractureProperty const& frac_prop = *_process_data.fracture_property;
    auto const& R = frac_prop.R;
    double const& dt = _process_data.dt;

    // the index of a normal (normal to a fracture plane) component
    // in a displacement vector
    auto const index_normal = GlobalDim - 1;

    typename ShapeMatricesTypePressure::NodalMatrixType laplace_p;
    laplace_p.setZero(pressure_size, pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage_p;
    storage_p.setZero(pressure_size, pressure_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        displacement_size, pressure_size>
        Kgp;
    Kgp.setZero(displacement_size, pressure_size);

    using GlobalDimMatrix = Eigen::Matrix<double, GlobalDim, GlobalDim>;
    using GlobalDimVector = Eigen::Matrix<double, GlobalDim, 1>;
    GlobalDimMatrix local_k_tensor;
    local_k_tensor.setZero();
    GlobalDimMatrix local_dk_db_tensor;
    local_dk_db_tensor.setZero();

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
        auto const& H_g = ip_data.H_u;
        auto const& identity2 =
            MaterialLib::Fracture::FractureIdentity2<GlobalDim>::value;

        auto& mat = ip_data.fracture_material;
        auto& effective_stress = ip_data.sigma_eff;
        auto const& effective_stress_prev = ip_data.sigma_eff_prev;
        auto& w = ip_data.w;
        auto const& w_prev = ip_data.w_prev;
        auto& C = ip_data.C;
        auto& b = ip_data.aperture;
        auto q = ip_data.darcy_velocity.head(GlobalDim);

        double const S = (*frac_prop.specific_storage)(t, x_position)[0];
        double const mu = _process_data.fluid_viscosity(t, x_position)[0];
        auto const alpha = (*frac_prop.biot_coefficient)(t, x_position)[0];
        auto const rho_fr = _process_data.fluid_density(t, x_position)[0];
        auto const& gravity_vec = _process_data.specific_body_force;

        // displacement jumps in local coordinates
        w.noalias() = R * H_g * nodal_g;

        // aperture
        b = ip_data.aperture0 + w[index_normal];

        // local C, local stress
        mat.computeConstitutiveRelation(
                    t, x_position,
                    w_prev, w,
                    effective_stress_prev, effective_stress, C);

        if (b < 1e-6) // < 0.0
        {
            //OGS_FATAL("Fracture aperture is %g, but it must be non-negative.", b);
            WARN("e %d, gp %d: Fracture aperture is %g, but it must be non-negative.", _element.getID(), ip, b);
            C(index_normal, index_normal) = 1e15;
            b = 1e-6;
        }

        if (_element.getID() == 55)
        {
            std::cout << "# e=" << _element.getID() << ", ip=" << ip << "\n";
            std::cout << "p=" << (N_p * nodal_p) << "\n";
            std::cout << "p_dot=" << (N_p * nodal_p_dot) << "\n";
            std::cout << "w=" << w.transpose() << "\n";
            std::cout << "sigma'_prev=" << effective_stress_prev.transpose() << "\n";
            std::cout << "sigma'=" << effective_stress.transpose() << "\n";
        }

        // permeability
        double const local_k = b * b / 12;
        ip_data.permeability = local_k;
        local_k_tensor.diagonal().head(GlobalDim-1).setConstant(local_k);
        GlobalDimMatrix const k = R.transpose() * local_k_tensor * R;

        // derivative of permeability respect to aperture
        double const local_dk_db = b / 6.;
        local_dk_db_tensor.diagonal().head(GlobalDim-1).setConstant(local_dk_db);
        GlobalDimMatrix const dk_db = R.transpose() * local_dk_db_tensor * R;

        // velocity
        GlobalDimVector const grad_head = dNdx_p * nodal_p + rho_fr * gravity_vec;
        q.noalias() = - k / mu * grad_head;

        //
        // displacement equation, displacement jump part
        //
        if (_process_data.use_initial_stress_as_reference)
        {
            auto const vec_sigma_eff_ref =
                _process_data.initial_fracture_effective_stress(t, x_position);
            auto const sigma_eff_ref =
                Eigen::Map<typename ShapeMatricesTypeDisplacement::
                               template VectorType<GlobalDim> const>(
                    vec_sigma_eff_ref.data(), GlobalDim);
            rhs_g.noalias() -= H_g.transpose() * R.transpose() *
                               (effective_stress - sigma_eff_ref) * ip_w;
        }
        else
        {
            rhs_g.noalias() -= H_g.transpose() * R.transpose() * effective_stress * ip_w;
        }
        J_gg.noalias() += H_g.transpose() * R.transpose() * C * R * H_g * ip_w;

        //
        // displacement equation, pressure part
        //
        Kgp.noalias() += H_g.transpose() * R.transpose() * alpha * identity2 * N_p * ip_w;

        //
        // pressure equation, pressure part.
        //
        storage_p.noalias() += N_p.transpose() * b * S * N_p * ip_w;
        laplace_p.noalias() += dNdx_p.transpose() * b * k / mu * dNdx_p * ip_w;
        rhs_p.noalias() += dNdx_p.transpose() * b * k / mu * rho_fr * gravity_vec * ip_w;

        //
        // pressure equation, displacement jump part.
        //
        Eigen::Matrix<double, 1, displacement_size> const mT_R_Hg = identity2.transpose() * R * H_g;
        J_pg.noalias() += N_p.transpose() * S * N_p * nodal_p_dot * mT_R_Hg * ip_w;
        J_pg.noalias() += - dNdx_p.transpose() * q * mT_R_Hg * ip_w;
        J_pg.noalias() += - dNdx_p.transpose() * b * (- dk_db / mu * grad_head) * mT_R_Hg * ip_w;
    }

    // displacement equation, pressure part
    J_gp.noalias() -= Kgp;

    // pressure equation, pressure part.
    J_pp.noalias() += laplace_p + storage_p / dt;

    // pressure equation, displacement jump part.
    J_pg.noalias() += Kgp.transpose() / dt;

    // pressure equation
    rhs_p.noalias() -=
        laplace_p * nodal_p + storage_p * nodal_p_dot + Kgp.transpose() * nodal_g_dot;

    // displacement equation
    if (_process_data.use_initial_stress_as_reference)
        rhs_g.noalias() -= - Kgp * (nodal_p - _initial_pressure);
    else
        rhs_g.noalias() -= - Kgp * nodal_p;
}


template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, unsigned GlobalDim>
void
HydroMechanicsLocalAssemblerFracture<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, IntegrationMethod,
                                     GlobalDim>::
computeSecondaryVariableConcreteWithVector(
                                const double t,
                                Eigen::VectorXd const& local_x)
{
    //auto const nodal_p = local_x.segment(pressure_index, pressure_size);
    auto const nodal_g = local_x.segment(displacement_index, displacement_size);

    FractureProperty const& frac_prop = *_process_data.fracture_property;
    auto const& R = frac_prop.R;
    // the index of a normal (normal to a fracture plane) component
    // in a displacement vector
    auto const index_normal = GlobalDim - 1;

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points = _ip_data.size();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto& ip_data = _ip_data[ip];
        auto const& H_g = ip_data.H_u;

        auto& mat = ip_data.fracture_material;
        auto& effective_stress = ip_data.sigma_eff;
        auto const& effective_stress_prev = ip_data.sigma_eff_prev;
        auto& w = ip_data.w;
        auto const& w_prev = ip_data.w_prev;
        auto& C = ip_data.C;
        auto& b = ip_data.aperture;

        // displacement jumps in local coordinates
        w.noalias() = R * H_g * nodal_g;

        // aperture
        b = ip_data.aperture0 + w[index_normal];

        // local C, local stress
        mat.computeConstitutiveRelation(
                    t, x_position,
                    w_prev, w,
                    effective_stress_prev, effective_stress, C);

        if (b < 1e-6) // < 0.0
        {
            //OGS_FATAL("Fracture aperture is %g, but it must be non-negative.", b);
            WARN("e %d, gp %d: Fracture aperture is %g, but it must be non-negative.", _element.getID(), ip, b);
        }

        // permeability
        double const local_k = b * b / 12;
        ip_data.permeability = local_k;
    }

    double ele_b = 0;
    double ele_k = 0;
    Eigen::Vector2d ele_w;
    ele_w.setZero();
    for (auto const& ip : _ip_data)
    {
        ele_b += ip.aperture;
        ele_k += ip.permeability;
        ele_w += ip.w;
    }
    ele_b /= _ip_data.size();
    ele_k /= _ip_data.size();
    ele_w /= _ip_data.size();
    (*_process_data.mesh_prop_b)[this->_element.getID()] = ele_b;
    (*_process_data.mesh_prop_k_f)[this->_element.getID()] = ele_k;
    (*_process_data.mesh_prop_w_n)[this->_element.getID()] = ele_w[index_normal];
    (*_process_data.mesh_prop_w_s)[this->_element.getID()] = ele_w[0];
}


template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, unsigned GlobalDim>
void
HydroMechanicsLocalAssemblerFracture<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, IntegrationMethod,
                                     GlobalDim>::
postTimestepConcrete(std::vector<double> const& /*local_x*/)
{
    double ele_b = 0;
    double ele_k = 0;
    for (auto const& ip : _ip_data)
    {
        ele_b += ip.aperture;
        ele_k += ip.permeability;
    }
    ele_b /= _ip_data.size();
    ele_k /= _ip_data.size();
    (*_process_data.mesh_prop_b)[this->_element.getID()] = ele_b;
    (*_process_data.mesh_prop_k_f)[this->_element.getID()] = ele_k;
}

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib

#endif  // PROCESSLIB_LIE_HYDROMECHANICS_HYDROMECHANICSLOCALASSEMBLER_FRACTURE_IMPL_H_
