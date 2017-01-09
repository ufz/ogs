/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "HydroMechanicsLocalAssemblerMatrixNearFracture.h"

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, unsigned GlobalDim>
HydroMechanicsLocalAssemblerMatrixNearFracture<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    GlobalDim>::
    HydroMechanicsLocalAssemblerMatrixNearFracture(
        MeshLib::Element const& e,
        std::size_t const n_variables,
        std::size_t const local_matrix_size,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        HydroMechanicsProcessData<GlobalDim>& process_data)
    : HydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                         ShapeFunctionPressure,
                                         IntegrationMethod, GlobalDim>(
          e, n_variables, local_matrix_size, dofIndex_to_localIndex,
          is_axially_symmetric, integration_order, process_data)
{
}


template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, unsigned GlobalDim>
void
HydroMechanicsLocalAssemblerMatrixNearFracture<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    GlobalDim>::
assembleWithJacobianConcrete(
    double const t,
    Eigen::VectorXd const& local_x,
    Eigen::VectorXd const& local_x_dot,
    Eigen::VectorXd& local_b,
    Eigen::MatrixXd& local_J)
{
    auto p = const_cast<Eigen::VectorXd&>(local_x).segment(pressure_index, pressure_size);
    auto p_dot = const_cast<Eigen::VectorXd&>(local_x_dot).segment(pressure_index, pressure_size);
    if (_process_data.deactivate_matrix_in_flow)
    {
        Base::setPressureOfInactiveNodes(t, p);
        Base::setPressureDotOfInactiveNodes(p_dot);
    }
    auto const u = local_x.segment(displacement_index, displacement_size);
    auto const u_dot =
        local_x_dot.segment(displacement_index, displacement_size);

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
    auto const& fracture_props = *_process_data.fracture_property;
    double const ele_levelset = calculateLevelSetFunction(
        fracture_props, _element.getCenterOfGravity().getCoords());

    if (ele_levelset == 0)
    {
        // no DoF exists for displacement jumps. do the normal assebmly
        Base::assembleBlockMatricesWithJacobian(t, p, p_dot, u, u_dot, rhs_p, rhs_u, J_pp, J_pu, J_uu, J_up);
        return;
    }

    // Displacement jumps should be taken into account

    // compute true displacements
    auto const g = local_x.segment(displacement_jump_index, displacement_size);
    auto const g_dot = local_x_dot.segment(displacement_jump_index, displacement_size);
    Eigen::VectorXd const total_u = u + ele_levelset * g;
    Eigen::VectorXd const total_u_dot = u_dot + ele_levelset * g_dot;

    // evaluate residuals and Jacobians for pressure and displacements
    Base::assembleBlockMatricesWithJacobian(t, p, p_dot, total_u, total_u_dot, rhs_p, rhs_u, J_pp, J_pu, J_uu, J_up);

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
          typename IntegrationMethod, unsigned GlobalDim>
void
HydroMechanicsLocalAssemblerMatrixNearFracture<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    GlobalDim>::
computeSecondaryVariableConcreteWithVector(
    double const t,
    Eigen::VectorXd const& local_x)
{
    auto p = const_cast<Eigen::VectorXd&>(local_x).segment(pressure_index, pressure_size);
    if (_process_data.deactivate_matrix_in_flow)
        Base::setPressureOfInactiveNodes(t, p);
    auto u = local_x.segment(displacement_index, displacement_size);

    // levelset value of the element
    // remark: this assumes the levelset function is uniform within an element
    auto const& fracture_props = *_process_data.fracture_property;
    double const ele_levelset = calculateLevelSetFunction(
        fracture_props, _element.getCenterOfGravity().getCoords());

    if (ele_levelset == 0)
    {
        // no DoF exists for displacement jumps. do the normal assebmly
        Base::computeSecondaryVariableConcreteWithBlockVectors(t, p, u);
        return;
    }

    // Displacement jumps should be taken into account

    // compute true displacements
    auto const g = local_x.segment(displacement_jump_index, displacement_size);
    Eigen::VectorXd const total_u = u + ele_levelset * g;

    // evaluate residuals and Jacobians for pressure and displacements
    Base::computeSecondaryVariableConcreteWithBlockVectors(t, p, total_u);
}


}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
