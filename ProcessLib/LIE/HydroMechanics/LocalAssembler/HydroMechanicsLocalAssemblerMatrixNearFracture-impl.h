/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "HydroMechanicsLocalAssemblerMatrixNearFracture.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int GlobalDim>
HydroMechanicsLocalAssemblerMatrixNearFracture<
    ShapeFunctionDisplacement, ShapeFunctionPressure, GlobalDim>::
    HydroMechanicsLocalAssemblerMatrixNearFracture(
        MeshLib::Element const& e,
        std::size_t const n_variables,
        std::size_t const local_matrix_size,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        HydroMechanicsProcessData<GlobalDim>& process_data)
    : HydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                         ShapeFunctionPressure, GlobalDim>(
          e, n_variables, local_matrix_size, dofIndex_to_localIndex,
          integration_method, is_axially_symmetric, process_data),
      _e_center_coords(getCenterOfGravity(e).data())
{
    // currently not supporting multiple fractures
    _fracture_props.push_back(process_data.fracture_property.get());
    _fracID_to_local.insert({0, 0});
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int GlobalDim>
void HydroMechanicsLocalAssemblerMatrixNearFracture<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    GlobalDim>::assembleWithJacobianConcrete(double const t, double const dt,
                                             Eigen::VectorXd const& local_x,
                                             Eigen::VectorXd const&
                                                 local_x_prev,
                                             Eigen::VectorXd& local_b,
                                             Eigen::MatrixXd& local_J)
{
    auto p = const_cast<Eigen::VectorXd&>(local_x).segment(pressure_index,
                                                           pressure_size);
    auto p_prev = const_cast<Eigen::VectorXd&>(local_x_prev)
                      .segment(pressure_index, pressure_size);
    if (_process_data.deactivate_matrix_in_flow)
    {
        Base::setPressureOfInactiveNodes(t, p);
    }
    auto const u = local_x.segment(displacement_index, displacement_size);
    auto const u_prev =
        local_x_prev.segment(displacement_index, displacement_size);

    auto rhs_p = local_b.segment(pressure_index, pressure_size);
    auto rhs_u = local_b.segment(displacement_index, displacement_size);

    auto J_pp = local_J.block(pressure_index, pressure_index, pressure_size,
                              pressure_size);
    auto J_pu = local_J.block(pressure_index, displacement_index, pressure_size,
                              displacement_size);
    auto J_up = local_J.block(displacement_index, pressure_index,
                              displacement_size, pressure_size);
    auto J_uu = local_J.block(displacement_index, displacement_index,
                              displacement_size, displacement_size);

    // levelset value of the element
    // remark: this assumes the levelset function is uniform within an element
    std::vector<double> levelsets = uGlobalEnrichments(
        _fracture_props, _junction_props, _fracID_to_local, _e_center_coords);
    double const ele_levelset = levelsets[0];  // single fracture

    if (ele_levelset == 0)
    {
        // no DoF exists for displacement jumps. do the normal assembly
        Base::assembleBlockMatricesWithJacobian(
            t, dt, p, p_prev, u, u_prev, rhs_p, rhs_u, J_pp, J_pu, J_uu, J_up);
        return;
    }

    // Displacement jumps should be taken into account

    // compute true displacements
    auto const g = local_x.segment(displacement_jump_index, displacement_size);
    auto const g_prev =
        local_x_prev.segment(displacement_jump_index, displacement_size);
    Eigen::VectorXd const total_u = u + ele_levelset * g;
    Eigen::VectorXd const total_u_prev = u_prev + ele_levelset * g_prev;

    // evaluate residuals and Jacobians for pressure and displacements
    Base::assembleBlockMatricesWithJacobian(t, dt, p, p_prev, total_u,
                                            total_u_prev, rhs_p, rhs_u, J_pp,
                                            J_pu, J_uu, J_up);

    // compute residuals and Jacobians for displacement jumps
    auto rhs_g = local_b.segment(displacement_jump_index, displacement_size);
    auto J_pg = local_J.block(pressure_index, displacement_jump_index,
                              pressure_size, displacement_size);
    auto J_ug = local_J.block(displacement_index, displacement_jump_index,
                              displacement_size, displacement_size);
    auto J_gp = local_J.block(displacement_jump_index, pressure_index,
                              displacement_size, pressure_size);
    auto J_gu = local_J.block(displacement_jump_index, displacement_index,
                              displacement_size, displacement_size);
    auto J_gg = local_J.block(displacement_jump_index, displacement_jump_index,
                              displacement_size, displacement_size);

    rhs_g = ele_levelset * rhs_u;
    J_pg = ele_levelset * J_pu;
    J_ug = ele_levelset * J_uu;
    J_gp = ele_levelset * J_up;
    J_gu = ele_levelset * J_uu;
    J_gg = ele_levelset * ele_levelset * J_uu;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int GlobalDim>
void HydroMechanicsLocalAssemblerMatrixNearFracture<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    GlobalDim>::postTimestepConcreteWithVector(double const t, double const dt,
                                               Eigen::VectorXd const& local_x)
{
    auto p = const_cast<Eigen::VectorXd&>(local_x).segment(pressure_index,
                                                           pressure_size);
    if (_process_data.deactivate_matrix_in_flow)
    {
        Base::setPressureOfInactiveNodes(t, p);
    }
    auto u = local_x.segment(displacement_index, displacement_size);

    // levelset value of the element
    // remark: this assumes the levelset function is uniform within an element
    std::vector<double> levelsets = uGlobalEnrichments(
        _fracture_props, _junction_props, _fracID_to_local, _e_center_coords);
    double const ele_levelset = levelsets[0];  // single fracture

    if (ele_levelset == 0)
    {
        // no DoF exists for displacement jumps. do the normal assembly
        Base::postTimestepConcreteWithBlockVectors(t, dt, p, u);
        return;
    }

    // Displacement jumps should be taken into account

    // compute true displacements
    auto const g = local_x.segment(displacement_jump_index, displacement_size);
    Eigen::VectorXd const total_u = u + ele_levelset * g;

    // evaluate residuals and Jacobians for pressure and displacements
    Base::postTimestepConcreteWithBlockVectors(t, dt, p, total_u);
}

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
