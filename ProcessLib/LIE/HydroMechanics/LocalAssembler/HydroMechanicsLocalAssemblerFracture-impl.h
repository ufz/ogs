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

#include <variant>

#include "HydroMechanicsLocalAssemblerFracture.h"
#include "MaterialLib/FractureModels/FractureIdentity2.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/Interpolation.h"
#include "ProcessLib/LIE/Common/LevelSetFunction.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{
namespace MPL = MaterialPropertyLib;

template <int DisplacementDim, typename RotationMatrix>
Eigen::Matrix<double, DisplacementDim, DisplacementDim> createRotatedTensor(
    RotationMatrix const& R, double const value)
{
    using M = Eigen::Matrix<double, DisplacementDim, DisplacementDim>;
    M tensor = M::Zero();
    tensor.diagonal().head(DisplacementDim - 1).setConstant(value);
    return (R.transpose() * tensor * R).eval();
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
HydroMechanicsLocalAssemblerFracture<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, DisplacementDim>::
    HydroMechanicsLocalAssemblerFracture(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        HydroMechanicsProcessData<DisplacementDim>& process_data)
    : HydroMechanicsLocalAssemblerInterface(
          e, is_axially_symmetric, integration_method,
          ShapeFunctionDisplacement::NPOINTS * DisplacementDim +
              ShapeFunctionPressure::NPOINTS,
          dofIndex_to_localIndex),
      _process_data(process_data)
{
    assert(e.getDimension() == DisplacementDim - 1);

    unsigned const n_integration_points =
        integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);
    _secondary_data.N.resize(n_integration_points);

    auto const shape_matrices_u =
        NumLib::initShapeMatrices<ShapeFunctionDisplacement,
                                  ShapeMatricesTypeDisplacement,
                                  DisplacementDim>(e, is_axially_symmetric,
                                                   integration_method);

    auto const shape_matrices_p =
        NumLib::initShapeMatrices<ShapeFunctionPressure,
                                  ShapeMatricesTypePressure, DisplacementDim>(
            e, is_axially_symmetric, integration_method);

    auto mat_id = (*_process_data.mesh_prop_materialIDs)[e.getID()];
    auto frac_id = _process_data.map_materialID_to_fractureID[mat_id];
    _fracture_property = &_process_data.fracture_properties[frac_id];

    // Get element nodes for aperture0 interpolation from nodes to integration
    // point. The aperture0 parameter is time-independent.
    typename ShapeMatricesTypeDisplacement::NodalVectorType
        aperture0_node_values =
            _fracture_property->aperture0.getNodalValuesOnElement(
                e, /*time independent*/ 0);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        _ip_data.emplace_back(*_process_data.fracture_model);
        auto const& sm_u = shape_matrices_u[ip];
        auto const& sm_p = shape_matrices_p[ip];
        ParameterLib::SpatialPosition x_position = {
            std::nullopt, _element.getID(),
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    _element, sm_u.N))};

        auto& ip_data = _ip_data[ip];
        ip_data.integration_weight =
            sm_u.detJ * sm_u.integralMeasure *
            integration_method.getWeightedPoint(ip).getWeight();

        ip_data.H_u.setZero(
            DisplacementDim,
            ShapeFunctionDisplacement::NPOINTS * DisplacementDim);
        computeHMatrix<
            DisplacementDim, ShapeFunctionDisplacement::NPOINTS,
            typename ShapeMatricesTypeDisplacement::NodalRowVectorType,
            HMatrixType>(sm_u.N, ip_data.H_u);
        ip_data.N_p = sm_p.N;
        ip_data.dNdx_p = sm_p.dNdx;

        _secondary_data.N[ip] = sm_u.N;

        // Initialize current time step values
        ip_data.w.setZero(DisplacementDim);
        ip_data.sigma_eff.setZero(DisplacementDim);

        // Previous time step values are not initialized and are set later.
        ip_data.w_prev.resize(DisplacementDim);
        ip_data.sigma_eff_prev.resize(DisplacementDim);

        ip_data.C.resize(DisplacementDim, DisplacementDim);

        ip_data.aperture0 = aperture0_node_values.dot(sm_u.N);
        ip_data.aperture = ip_data.aperture0;

        auto const initial_effective_stress =
            _process_data.initial_fracture_effective_stress(0, x_position);
        for (int i = 0; i < DisplacementDim; i++)
        {
            ip_data.sigma_eff[i] = initial_effective_stress[i];
            ip_data.sigma_eff_prev[i] = initial_effective_stress[i];
        }
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void HydroMechanicsLocalAssemblerFracture<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    assembleWithJacobianConcrete(double const t, double const dt,
                                 Eigen::VectorXd const& local_x,
                                 Eigen::VectorXd const& local_x_prev,
                                 Eigen::VectorXd& local_b,
                                 Eigen::MatrixXd& local_J)
{
    auto const p = local_x.segment(pressure_index, pressure_size);
    auto const p_prev = local_x_prev.segment(pressure_index, pressure_size);
    auto const g = local_x.segment(displacement_index, displacement_size);
    auto const g_prev =
        local_x_prev.segment(displacement_index, displacement_size);

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

    assembleBlockMatricesWithJacobian(t, dt, p, p_prev, g, g_prev, rhs_p, rhs_g,
                                      J_pp, J_pg, J_gg, J_gp);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void HydroMechanicsLocalAssemblerFracture<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    assembleBlockMatricesWithJacobian(
        double const t, double const dt,
        Eigen::Ref<const Eigen::VectorXd> const& p,
        Eigen::Ref<const Eigen::VectorXd> const& p_prev,
        Eigen::Ref<const Eigen::VectorXd> const& g,
        Eigen::Ref<const Eigen::VectorXd> const& g_prev,
        Eigen::Ref<Eigen::VectorXd> rhs_p, Eigen::Ref<Eigen::VectorXd> rhs_g,
        Eigen::Ref<Eigen::MatrixXd> J_pp, Eigen::Ref<Eigen::MatrixXd> J_pg,
        Eigen::Ref<Eigen::MatrixXd> J_gg, Eigen::Ref<Eigen::MatrixXd> J_gp)
{
    auto const& R = _fracture_property->R;

    // the index of a normal (normal to a fracture plane) component
    // in a displacement vector
    auto constexpr index_normal = DisplacementDim - 1;

    typename ShapeMatricesTypePressure::NodalMatrixType laplace_p =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage_p =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        displacement_size, pressure_size>
        Kgp = ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size, pressure_size>::Zero(displacement_size,
                                                    pressure_size);

    // Projection of the body force vector at the element.
    Eigen::MatrixXd const global2local_rotation =
        R.template topLeftCorner<ShapeFunctionPressure::DIM, DisplacementDim>();

    GlobalDimVectorType const gravity_vec = global2local_rotation.transpose() *
                                            global2local_rotation *
                                            _process_data.specific_body_force;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    MPL::VariableArray variables;
    auto const& medium = _process_data.media_map.getMedium(_element.getID());
    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const T_ref =
        medium->property(MPL::PropertyType::reference_temperature)
            .template value<double>(variables, x_position, t, dt);
    variables.temperature = T_ref;

    unsigned const n_integration_points = _ip_data.size();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto& ip_data = _ip_data[ip];
        auto const& ip_w = ip_data.integration_weight;
        auto const& N_p = ip_data.N_p;
        auto const& dNdx_p = ip_data.dNdx_p;
        auto const& H_g = ip_data.H_u;
        auto const& identity2 =
            MaterialLib::Fracture::FractureIdentity2<DisplacementDim>::value;

        x_position = {
            std::nullopt, _element.getID(),
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionPressure,
                                               ShapeMatricesTypePressure>(
                    _element, N_p))};

        auto& mat = ip_data.fracture_material;
        auto& effective_stress = ip_data.sigma_eff;
        auto const& effective_stress_prev = ip_data.sigma_eff_prev;
        auto& w = ip_data.w;
        auto const& w_prev = ip_data.w_prev;
        auto& C = ip_data.C;
        auto& state = *ip_data.material_state_variables;
        auto& b_m = ip_data.aperture;

        auto const rho_fr =
            liquid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);
        variables.density = rho_fr;

        auto const alpha =
            medium->property(MPL::PropertyType::biot_coefficient)
                .template value<double>(variables, x_position, t, dt);

        double const S =
            medium->property(MPL::PropertyType::storage)
                .template value<double>(variables, x_position, t, dt);

        auto const mu =
            liquid_phase.property(MPL::PropertyType::viscosity)
                .template value<double>(variables, x_position, t, dt);

        // displacement jumps in local coordinates
        w.noalias() = R * H_g * g;

        // aperture
        b_m = ip_data.aperture0 + w[index_normal];
        if (b_m < 0.0)
        {
            DBUG(
                "Element {:d}, gp {:d}: Fracture aperture is {:g}, but it must "
                "be "
                "non-negative. Setting it to zero.",
                _element.getID(), ip, b_m);
            b_m = 0;
        }

        auto const initial_effective_stress =
            _process_data.initial_fracture_effective_stress(0, x_position);

        Eigen::Map<typename HMatricesType::ForceVectorType const> const stress0(
            initial_effective_stress.data(), initial_effective_stress.size());

        // local C, local stress
        mat.computeConstitutiveRelation(
            t, x_position, ip_data.aperture0, stress0, w_prev, w,
            effective_stress_prev, effective_stress, C, state);

        //
        // displacement equation, displacement jump part
        //
        rhs_g.noalias() -=
            H_g.transpose() * R.transpose() * effective_stress * ip_w;
        J_gg.noalias() += H_g.transpose() * R.transpose() * C * R * H_g * ip_w;

        //
        // displacement equation, pressure part
        //
        Kgp.noalias() +=
            H_g.transpose() * R.transpose() * alpha * identity2 * N_p * ip_w;

        //
        // pressure equation, pressure part.
        //

        variables.fracture_aperture = b_m;
        // Assume that the fracture permeability is isotropic
        auto const permeability =
            medium->property(MPL::PropertyType::permeability)
                .value(variables, x_position, t, dt);

        auto& k = ip_data.permeability;
        k = std::get<double>(permeability);
        double const k_over_mu = k / mu;
        storage_p.noalias() += N_p.transpose() * b_m * S * N_p * ip_w;
        laplace_p.noalias() +=
            dNdx_p.transpose() * b_m * k_over_mu * dNdx_p * ip_w;
        rhs_p.noalias() +=
            dNdx_p.transpose() * b_m * k_over_mu * rho_fr * gravity_vec * ip_w;

        //
        // pressure equation, displacement jump part.
        //
        GlobalDimVectorType const grad_head = dNdx_p * p + rho_fr * gravity_vec;
        Eigen::Matrix<double, 1, displacement_size> const mT_R_Hg =
            identity2.transpose() * R * H_g;
        // velocity in global coordinates
        ip_data.darcy_velocity = -k_over_mu * grad_head;
        J_pg.noalias() +=
            N_p.transpose() * S * N_p * (p - p_prev) / dt * mT_R_Hg * ip_w;

        // derivative of permeability with respect to aperture
        double const dk_db_over_mu =
            medium->property(MPL::PropertyType::permeability)
                .template dValue<double>(variables,
                                         MPL::Variable::fracture_aperture,
                                         x_position, t, dt) /
            mu;
        J_pg.noalias() +=
            dNdx_p.transpose() * k_over_mu * grad_head * mT_R_Hg * ip_w;
        J_pg.noalias() += dNdx_p.transpose() * b_m * dk_db_over_mu * grad_head *
                          mT_R_Hg * ip_w;
    }

    // displacement equation, pressure part
    J_gp.noalias() -= Kgp;

    // pressure equation, pressure part.
    J_pp.noalias() += laplace_p + storage_p / dt;

    // pressure equation, displacement jump part.
    J_pg.noalias() += Kgp.transpose() / dt;

    // pressure equation
    rhs_p.noalias() -= laplace_p * p + storage_p * (p - p_prev) / dt +
                       Kgp.transpose() * (g - g_prev) / dt;

    // displacement equation
    rhs_g.noalias() -= -Kgp * p;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void HydroMechanicsLocalAssemblerFracture<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::postTimestepConcreteWithVector(const double t,
                                                     double const /*dt*/,
                                                     Eigen::VectorXd const&
                                                         local_x)
{
    auto const nodal_g = local_x.segment(displacement_index, displacement_size);

    auto const& R = _fracture_property->R;
    // the index of a normal (normal to a fracture plane) component
    // in a displacement vector
    auto constexpr index_normal = DisplacementDim - 1;

    ParameterLib::SpatialPosition x_position;
    auto const e_id = _element.getID();
    x_position.setElementID(e_id);

    unsigned const n_integration_points = _ip_data.size();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto& ip_data = _ip_data[ip];
        auto const& H_g = ip_data.H_u;
        auto& mat = ip_data.fracture_material;
        auto& effective_stress = ip_data.sigma_eff;
        auto const& effective_stress_prev = ip_data.sigma_eff_prev;
        auto& w = ip_data.w;
        auto const& w_prev = ip_data.w_prev;
        auto& C = ip_data.C;
        auto& state = *ip_data.material_state_variables;
        auto& b_m = ip_data.aperture;

        x_position = {
            std::nullopt, e_id,
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionPressure,
                                               ShapeMatricesTypePressure>(
                    _element, ip_data.N_p))};

        // displacement jumps in local coordinates
        w.noalias() = R * H_g * nodal_g;

        // aperture
        b_m = ip_data.aperture0 + w[index_normal];
        if (b_m < 0.0)
        {
            DBUG(
                "Element {:d}, gp {:d}: Fracture aperture is {:g}, but it is "
                "expected to be non-negative. Setting it to zero now.",
                _element.getID(), ip, b_m);
            b_m = 0;
        }

        auto const initial_effective_stress =
            _process_data.initial_fracture_effective_stress(0, x_position);

        Eigen::Map<typename HMatricesType::ForceVectorType const> const stress0(
            initial_effective_stress.data(), initial_effective_stress.size());

        // local C, local stress
        mat.computeConstitutiveRelation(
            t, x_position, ip_data.aperture0, stress0, w_prev, w,
            effective_stress_prev, effective_stress, C, state);
    }

    double ele_b = 0;
    double ele_k = 0;
    typename HMatricesType::ForceVectorType ele_sigma_eff =
        HMatricesType::ForceVectorType::Zero(DisplacementDim);
    typename HMatricesType::ForceVectorType ele_w =
        HMatricesType::ForceVectorType::Zero(DisplacementDim);
    GlobalDimVectorType ele_velocity = GlobalDimVectorType::Zero();

    double ele_Fs = -std::numeric_limits<double>::max();
    for (auto const& ip : _ip_data)
    {
        ele_b += ip.aperture;
        ele_k += ip.permeability;
        ele_w += ip.w;
        ele_sigma_eff += ip.sigma_eff;
        ele_velocity += ip.darcy_velocity;
        ele_Fs = std::max(
            ele_Fs, ip.material_state_variables->getShearYieldFunctionValue());
    }
    ele_b /= static_cast<double>(n_integration_points);
    ele_k /= static_cast<double>(n_integration_points);
    ele_w /= static_cast<double>(n_integration_points);
    ele_sigma_eff /= static_cast<double>(n_integration_points);
    ele_velocity /= static_cast<double>(n_integration_points);
    auto const element_id = _element.getID();
    (*_process_data.mesh_prop_b)[element_id] = ele_b;
    (*_process_data.mesh_prop_k_f)[element_id] = ele_k;

    Eigen::Map<GlobalDimVectorType>(
        &(*_process_data.element_fracture_stresses)[e_id * DisplacementDim]) =
        ele_sigma_eff;

    Eigen::Map<GlobalDimVectorType>(
        &(*_process_data.element_fracture_velocities)[e_id * DisplacementDim]) =
        ele_velocity;

    Eigen::Map<GlobalDimVectorType>(
        &(*_process_data.element_local_jumps)[e_id * DisplacementDim]) = ele_w;

    (*_process_data.mesh_prop_fracture_shear_failure)[element_id] = ele_Fs;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> const& HydroMechanicsLocalAssemblerFracture<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    getIntPtFractureVelocity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    unsigned const n_integration_points = _ip_data.size();
    cache.clear();
    auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, DisplacementDim, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        cache_matrix.col(ip).noalias() = _ip_data[ip].darcy_velocity;
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> const& HydroMechanicsLocalAssemblerFracture<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    getIntPtFractureStress(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    unsigned const n_integration_points = _ip_data.size();
    cache.clear();
    auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, DisplacementDim, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        cache_matrix.col(ip).noalias() = _ip_data[ip].sigma_eff;
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> const& HydroMechanicsLocalAssemblerFracture<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    getIntPtFractureAperture(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(
        _ip_data, &IntegrationPointDataType::aperture, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> const& HydroMechanicsLocalAssemblerFracture<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    getIntPtFracturePermeability(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(
        _ip_data, &IntegrationPointDataType::permeability, cache);
}

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
