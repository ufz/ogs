/**
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on January 8, 2018, 3:00 PM
 */
#pragma once

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <utility>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/MPL/Utils/GetSymmetricTensor.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
namespace HMPhaseField
{
template <typename ShapeFunction, int DisplacementDim>
void HMPhaseFieldLocalAssembler<ShapeFunction, DisplacementDim>::
    assembleWithJacobianForStaggeredScheme(double const t, double const dt,
                                           Eigen::VectorXd const& local_x,
                                           Eigen::VectorXd const& local_x_prev,
                                           int const process_id,
                                           std::vector<double>& local_b_data,
                                           std::vector<double>& local_Jac_data)
{
    // For the equations with phase field.
    if (process_id == _process_data._phasefield_process_id)
    {
        assembleWithJacobianPhaseFieldEquations(t, dt, local_x, local_b_data,
                                                local_Jac_data);
        return;
    }

    // For the equations for hydro
    if (process_id == _process_data._hydro_process_id)
    {
        assembleWithJacobianHydroEquations(t, dt, local_x, local_x_prev,
                                           local_b_data, local_Jac_data);
        return;
    }

    // For the equations with deformation
    assembleWithJacobianForDeformationEquations(t, dt, local_x, local_b_data,
                                                local_Jac_data);
}

template <typename ShapeFunction, int DisplacementDim>
void HMPhaseFieldLocalAssembler<ShapeFunction, DisplacementDim>::
    assembleWithJacobianForDeformationEquations(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    auto const d = local_x.template segment<phasefield_size>(phasefield_index);
    auto const u =
        local_x.template segment<displacement_size>(displacement_index);
    auto const p = local_x.template segment<pressure_size>(pressure_index);

    auto local_Jac = MathLib::createZeroedMatrix<DeformationMatrix>(
        local_Jac_data, displacement_size, displacement_size);

    auto local_rhs = MathLib::createZeroedVector<DeformationVector>(
        local_b_data, displacement_size);

    auto const& medium = _process_data.media_map.getMedium(_element.getID());
    auto const& solid = medium->phase("Solid");
    auto const& fluid = fluidPhase(*medium);
    MPL::VariableArray vars;

    auto const& identity2 = Invariants::identity2;

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;
        double const d_ip = N.dot(d);
        double const p_ip = N.dot(p);

        ParameterLib::SpatialPosition const x_position{
            std::nullopt, this->_element.getID(),
            MathLib::Point3d(NumLib::interpolateCoordinates<ShapeFunction,
                                                            ShapeMatricesType>(
                this->_element, N))};

        auto const x_coord =
            x_position.getCoordinates().value()[0];  // r for axisymmetry

        auto const& B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunction::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx, N, x_coord, _is_axially_symmetric);

        auto& eps = _ip_data[ip].eps;
        eps.noalias() = B * u;

        double const k = _process_data.residual_stiffness(t, x_position)[0];
        double const ls = _process_data.crack_length_scale(t, x_position)[0];

        double const degradation =
            _process_data.degradation_derivative->degradation(d_ip, k, ls);
        _ip_data[ip].updateConstitutiveRelation(
            t, x_position, dt, u, degradation,
            _process_data.energy_split_model);

        auto const& sigma = _ip_data[ip].sigma;
        auto const& D = _ip_data[ip].D;

        auto& biot_coefficient = _ip_data[ip].biot_coefficient;
        auto& biot_modulus_inv = _ip_data[ip].biot_modulus_inv;
        auto const& fracture_enhanced_porosity =
            _ip_data[ip].fracture_enhanced_porosity;

        // Update the effective bulk modulus
        auto const& P_sph = Invariants::spherical_projection;
        auto const D_sph = P_sph * D * identity2;
        double const bulk_modulus_eff = Invariants::trace(D_sph) / 9.;

        auto const& solid_material =
            MaterialLib::Solids::selectSolidConstitutiveRelation(
                _process_data.solid_materials, _process_data.material_ids,
                _element.getID());
        auto const bulk_modulus = solid_material.getBulkModulus(t, x_position);
        auto const alpha_0 =
            solid.property(MPL::PropertyType::biot_coefficient)
                .template value<double>(vars, x_position, t, dt);
        auto const porosity_0 =
            medium->property(MPL::PropertyType::porosity)
                .template value<double>(vars, x_position, t, dt);
        double const bulk_modulus_degradation = bulk_modulus_eff / bulk_modulus;

        // Update Biot's coefficient
        biot_coefficient = 1. - bulk_modulus_degradation * (1. - alpha_0);

        // The reference porosity
        auto const porosity_reference = porosity_0 + fracture_enhanced_porosity;

        // Update Biot's modulus
        biot_modulus_inv = bulk_modulus_degradation * (alpha_0 - porosity_0) *
                           (1. - alpha_0) / bulk_modulus;

        auto const rho_sr =
            solid.property(MPL::PropertyType::density)
                .template value<double>(vars, x_position, t, dt);
        auto const rho_fr =
            fluid.property(MPL::PropertyType::density)
                .template value<double>(vars, x_position, t, dt);

        auto const rho =
            rho_sr * (1. - porosity_reference) + porosity_reference * rho_fr;
        auto const& b = _process_data.specific_body_force;

        local_rhs.noalias() -=
            (B.transpose() * (sigma - biot_coefficient * p_ip * identity2) -
             N_u_op(N).transpose() * rho * b) *
            w;
        local_Jac.noalias() += B.transpose() * D * B * w;
    }
}

template <typename ShapeFunction, int DisplacementDim>
void HMPhaseFieldLocalAssembler<ShapeFunction, DisplacementDim>::
    assembleWithJacobianHydroEquations(double const t, double const dt,
                                       Eigen::VectorXd const& local_x,
                                       Eigen::VectorXd const& local_x_prev,
                                       std::vector<double>& local_b_data,
                                       std::vector<double>& local_Jac_data)
{
    auto const d = local_x.template segment<phasefield_size>(phasefield_index);

    auto const p = local_x.template segment<pressure_size>(pressure_index);
    auto const p_prev =
        local_x_prev.template segment<pressure_size>(pressure_index);

    auto local_Jac = MathLib::createZeroedMatrix<PressureMatrix>(
        local_Jac_data, pressure_size, pressure_size);

    auto local_rhs = MathLib::createZeroedVector<PressureVector>(local_b_data,
                                                                 pressure_size);

    typename ShapeMatricesType::NodalMatrixType mass =
        ShapeMatricesType::NodalMatrixType::Zero(pressure_size, pressure_size);

    typename ShapeMatricesType::NodalMatrixType laplace =
        ShapeMatricesType::NodalMatrixType::Zero(pressure_size, pressure_size);

    typename ShapeMatricesType::NodalMatrixType stablizing =
        ShapeMatricesType::NodalMatrixType::Zero(pressure_size, pressure_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    auto const& medium = _process_data.media_map.getMedium(_element.getID());
    auto const& fluid = fluidPhase(*medium);
    MPL::VariableArray vars;

    auto const& P_sph = Invariants::spherical_projection;
    auto const& P_dev = Invariants::deviatoric_projection;
    auto const& identity2 = Invariants::identity2;
    auto const& ones2 = Invariants::ones2;

    double const k = _process_data.residual_stiffness(t, x_position)[0];
    double const ls = _process_data.crack_length_scale(t, x_position)[0];
    double const width = (*_process_data.width)[_element.getID()];
    double const fracture_threshold = _process_data.fracture_threshold;
    double const fracture_permeability_parameter =
        _process_data.fracture_permeability_parameter;
    double const fixed_stress_stabilization_parameter =
        _process_data.fixed_stress_stabilization_parameter;
    double const spatial_stabilization_parameter =
        _process_data.spatial_stabilization_parameter;
    auto const he =
        ls / _process_data.diffused_range_parameter;  // element size

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;
        double const d_ip = N.dot(d);

        auto const rho_fr =
            fluid.property(MPL::PropertyType::density)
                .template value<double>(vars, x_position, t, dt);
        double const cf = _process_data.fluid_compressibility;
        auto const Km = medium->property(MPL::PropertyType::permeability)
                            .template value<double>(vars, x_position, t, dt);
        auto const K = MPL::formEigenTensor<DisplacementDim>(Km);
        auto const mu = fluid.property(MPL::PropertyType::viscosity)
                            .template value<double>(vars, x_position, t, dt);

        auto const vol_strain = Invariants::trace(_ip_data[ip].eps);
        auto const vol_strain_prev = Invariants::trace(_ip_data[ip].eps_prev);
        auto const& biot_coefficient = _ip_data[ip].biot_coefficient;
        auto const& biot_modulus_inv = _ip_data[ip].biot_modulus_inv;
        auto const& fracture_enhanced_porosity =
            _ip_data[ip].fracture_enhanced_porosity;

        // The reference porosity
        auto const porosity_0 =
            medium->property(MPL::PropertyType::porosity)
                .template value<double>(vars, x_position, t, dt);
        auto const porosity_reference = porosity_0 + fracture_enhanced_porosity;

        double const dv_dt = (vol_strain - vol_strain_prev) / dt;

        // The degraded shear modulus
        auto const& D = _ip_data[ip].D;
        auto const D_sph = P_sph * D * identity2;
        auto const D_dev = P_dev * D * (ones2 - identity2) / std::sqrt(2.);
        auto const degraded_shear_modulus =
            Invariants::FrobeniusNorm(D_dev) / 2.;
        auto const degraded_bulk_modulus = Invariants::trace(D_sph) / 9.;

        double const residual_bulk_modulus = [&]
        {
            if ((*_process_data.ele_d)[_element.getID()] <
                _process_data.fracture_threshold)
            {
                // The residual bulk modulus in the fractured element
                double const degradation_threshold =
                    _process_data.degradation_derivative->degradation(
                        fracture_threshold, k, ls);
                auto const& D_threshold =
                    degradation_threshold * _ip_data[ip].C_tensile +
                    _ip_data[ip].C_compressive;
                auto const D_sph_threshold = P_sph * D_threshold * identity2;

                return Invariants::trace(D_sph_threshold) / 9.;
            }
            return degraded_bulk_modulus;
        }();

        double const modulus_rm = fixed_stress_stabilization_parameter *
                                  biot_coefficient * biot_coefficient /
                                  residual_bulk_modulus;

        double const stablization_spatial =
            spatial_stabilization_parameter * 0.25 * he * he /
            (degraded_bulk_modulus + 4. / 3. * degraded_shear_modulus);
        stablizing.noalias() +=
            dNdx.transpose() * stablization_spatial * dNdx * w;

        mass.noalias() +=
            (biot_modulus_inv + cf * porosity_reference + modulus_rm) *
            N.transpose() * N * w;

        auto const K_over_mu = K / mu;
        laplace.noalias() += dNdx.transpose() * K_over_mu * dNdx * w;
        auto const& b = _process_data.specific_body_force;

        // bodyforce-driven Darcy flow
        local_rhs.noalias() += dNdx.transpose() * rho_fr * K_over_mu * b * w;

        local_rhs.noalias() -= (biot_coefficient * dv_dt -
                                modulus_rm * _ip_data[ip].coupling_pressure) *
                               N.transpose() * w;

        if ((*_process_data.ele_d)[_element.getID()] <
            _process_data.irreversible_threshold)
        {
            // Fracture-enhanced permeability
            auto const& normal_ip = _ip_data[ip].normal_ip;
            auto const Kf =
                std::pow(1. - d_ip, fracture_permeability_parameter) * width *
                width * width / he / 12.0 *
                (Eigen::Matrix<double, DisplacementDim,
                               DisplacementDim>::Identity() -
                 normal_ip * normal_ip.transpose());
            laplace.noalias() += dNdx.transpose() * Kf / mu * dNdx * w;
        }
    }
    local_Jac.noalias() = laplace + mass / dt + stablizing / dt;

    local_rhs.noalias() -=
        laplace * p + mass * (p - p_prev) / dt + stablizing * (p - p_prev) / dt;
}

template <typename ShapeFunction, int DisplacementDim>
void HMPhaseFieldLocalAssembler<ShapeFunction, DisplacementDim>::
    assembleWithJacobianPhaseFieldEquations(double const t, double const dt,
                                            Eigen::VectorXd const& local_x,
                                            std::vector<double>& local_b_data,
                                            std::vector<double>& local_Jac_data)
{
    auto const d = local_x.template segment<phasefield_size>(phasefield_index);
    auto const p = local_x.template segment<pressure_size>(pressure_index);
    auto const u =
        local_x.template segment<displacement_size>(displacement_index);

    auto local_Jac = MathLib::createZeroedMatrix<PhaseFieldMatrix>(
        local_Jac_data, phasefield_size, phasefield_size);
    auto local_rhs = MathLib::createZeroedVector<PhaseFieldVector>(
        local_b_data, phasefield_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    auto const& solid_material =
        MaterialLib::Solids::selectSolidConstitutiveRelation(
            _process_data.solid_materials, _process_data.material_ids,
            _element.getID());

    auto const bulk_modulus = solid_material.getBulkModulus(t, x_position);
    auto const& medium = _process_data.media_map.getMedium(_element.getID());
    auto const& solid = medium->phase("Solid");
    MPL::VariableArray vars;

    double const k = _process_data.residual_stiffness(t, x_position)[0];
    double const ls = _process_data.crack_length_scale(t, x_position)[0];
    double const gc = _process_data.crack_resistance(t, x_position)[0];

    auto const& identity2 = Invariants::identity2;

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        double const d_ip = N.dot(d);
        double const p_ip = N.dot(p);
        double const degradation =
            _process_data.degradation_derivative->degradation(d_ip, k, ls);
        double const degradation_df1 =
            _process_data.degradation_derivative->degradationDf1(d_ip, k, ls);
        double const degradation_df2 =
            _process_data.degradation_derivative->degradationDf2(d_ip, k, ls);

        _ip_data[ip].updateConstitutiveRelation(
            t, x_position, dt, u, degradation,
            _process_data.energy_split_model);

        auto& biot_coefficient = _ip_data[ip].biot_coefficient;
        auto& biot_modulus_inv = _ip_data[ip].biot_modulus_inv;

        auto const alpha_0 =
            solid.property(MPL::PropertyType::biot_coefficient)
                .template value<double>(vars, x_position, t, dt);
        auto const porosity_0 =
            medium->property(MPL::PropertyType::porosity)
                .template value<double>(vars, x_position, t, dt);

        auto const& D = _ip_data[ip].D;
        auto const& P_sph = Invariants::spherical_projection;
        auto const D_sph = P_sph * D * identity2;
        double const bulk_modulus_eff = Invariants::trace(D_sph) / 9.;
        double const bulk_modulus_degradation = bulk_modulus_eff / bulk_modulus;

        // Update Biot's coefficient
        biot_coefficient = 1. - bulk_modulus_degradation * (1. - alpha_0);

        // Update Biot's modulus
        biot_modulus_inv = bulk_modulus_degradation * (alpha_0 - porosity_0) *
                           (1. - alpha_0) / bulk_modulus;

        auto const& strain_energy_tensile = _ip_data[ip].strain_energy_tensile;
        auto const& C_tensile = _ip_data[ip].C_tensile;
        auto const C_tensile_sph = P_sph * C_tensile * identity2;
        double const bulk_modulus_plus = Invariants::trace(C_tensile_sph) / 9.;

        auto const driven_energy =
            N.transpose() *
            (strain_energy_tensile + p_ip * p_ip / 2. * bulk_modulus_plus /
                                         bulk_modulus * (alpha_0 - porosity_0) *
                                         (1. - alpha_0) / bulk_modulus) *
            w;

        local_Jac.noalias() += driven_energy * N * degradation_df2;

        local_rhs.noalias() -= driven_energy * degradation_df1;

        calculateCrackLocalJacobianAndResidual<
            decltype(dNdx), decltype(N), decltype(w), decltype(d),
            decltype(local_Jac), decltype(local_rhs)>(
            dNdx, N, w, d, local_Jac, local_rhs, gc, ls,
            _process_data.phasefield_model);
    }
}

template <typename ShapeFunction, int DisplacementDim>
void HMPhaseFieldLocalAssembler<ShapeFunction, DisplacementDim>::
    postNonLinearSolverConcrete(Eigen::VectorXd const& local_x,
                                Eigen::VectorXd const& local_x_prev,
                                double const /*t*/, double const dt,
                                int const /*process_id*/)
{
    int const n_integration_points = _integration_method.getNumberOfPoints();
    auto const p = local_x.template segment<pressure_size>(pressure_index);
    auto const p_prev =
        local_x_prev.template segment<pressure_size>(pressure_index);

    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const& N = _ip_data[ip].N;
        _ip_data[ip].coupling_pressure = N.dot(p - p_prev) / dt;
    }
}

template <typename ShapeFunction, int DisplacementDim>
void HMPhaseFieldLocalAssembler<ShapeFunction, DisplacementDim>::
    approximateFractureWidth(
        std::size_t mesh_item_id,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
        std::vector<GlobalVector*> const& x, double const t,
        double const /*dt*/)
{
    std::vector<std::vector<GlobalIndexType>> indices_of_processes;
    indices_of_processes.reserve(dof_tables.size());
    std::transform(dof_tables.begin(), dof_tables.end(),
                   std::back_inserter(indices_of_processes),
                   [&](auto const dof_table)
                   { return NumLib::getIndices(mesh_item_id, *dof_table); });

    auto local_coupled_xs = getCoupledLocalSolutions(x, indices_of_processes);
    assert(local_coupled_xs.size() ==
           phasefield_size + displacement_size + pressure_size);

    auto const d = Eigen::Map<PhaseFieldVector const>(
        &local_coupled_xs[phasefield_index], phasefield_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    double const ele_d = std::clamp(d.sum() / d.size(), 0.0, 1.0);
    (*_process_data.ele_d)[_element.getID()] = ele_d;

    if ((*_process_data.ele_d)[_element.getID()] <
        _process_data.irreversible_threshold)
    {
        double const width_init = _process_data.width_init(t, x_position)[0];
        double const k = _process_data.residual_stiffness(t, x_position)[0];
        double const ls = _process_data.crack_length_scale(t, x_position)[0];
        double const he = ls / _process_data.diffused_range_parameter;

        int const n_integration_points =
            _integration_method.getNumberOfPoints();
        double width = 0.0;
        for (int ip = 0; ip < n_integration_points; ip++)
        {
            auto eps_tensor =
                MathLib::KelvinVector::kelvinVectorToTensor(_ip_data[ip].eps);
            Eigen::EigenSolver<decltype(eps_tensor)> eigen_solver(eps_tensor);
            Eigen::MatrixXf::Index maxIndex;
            double const max_principal_strain =
                eigen_solver.eigenvalues().real().maxCoeff(&maxIndex);
            auto const max_eigen_vector =
                eigen_solver.eigenvectors().real().col(maxIndex);

            // Fracture aperture estimation
            auto& width_ip = _ip_data[ip].width_ip;
            width_ip = max_principal_strain * he;
            width_ip = width_ip < width_init ? width_init : width_ip;
            width += width_ip;

            // Fracture direction estimation
            auto& normal_ip = _ip_data[ip].normal_ip;
            if (std::abs(max_principal_strain) > k)
            {
                for (int i = 0; i < DisplacementDim; i++)
                {
                    normal_ip[i] = max_eigen_vector[i];
                }
            }

            // Fracture enhanced porosity
            auto& fracture_enhanced_porosity =
                _ip_data[ip].fracture_enhanced_porosity;
            fracture_enhanced_porosity = width_ip / he;
        }

        // Update aperture for the fractured element
        (*_process_data.width)[_element.getID()] = width / n_integration_points;
    }
}

template <typename ShapeFunction, int DisplacementDim>
void HMPhaseFieldLocalAssembler<ShapeFunction, DisplacementDim>::computeEnergy(
    std::size_t mesh_item_id,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
    std::vector<GlobalVector*> const& x, double const t, double& elastic_energy,
    double& surface_energy, double& pressure_work)
{
    std::vector<std::vector<GlobalIndexType>> indices_of_processes;
    indices_of_processes.reserve(dof_tables.size());
    std::transform(dof_tables.begin(), dof_tables.end(),
                   std::back_inserter(indices_of_processes),
                   [&](auto const dof_table)
                   { return NumLib::getIndices(mesh_item_id, *dof_table); });

    auto const local_coupled_xs =
        getCoupledLocalSolutions(x, indices_of_processes);
    assert(local_coupled_xs.size() ==
           phasefield_size + displacement_size + pressure_size);

    auto const d = Eigen::Map<PhaseFieldVector const>(
        &local_coupled_xs[phasefield_index], phasefield_size);
    auto const u = Eigen::Map<DeformationVector const>(
        &local_coupled_xs[displacement_index], displacement_size);
    auto const p = Eigen::Map<PressureVector const>(
        &local_coupled_xs[pressure_index], pressure_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    double element_elastic_energy = 0.0;
    double element_surface_energy = 0.0;
    double element_pressure_work = 0.0;

    double const gc = _process_data.crack_resistance(t, x_position)[0];
    double const ls = _process_data.crack_length_scale(t, x_position)[0];
    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;
        double const d_ip = N.dot(d);
        double const p_ip = N.dot(p);

        element_elastic_energy += _ip_data[ip].elastic_energy * w;

        switch (_process_data.phasefield_model)
        {
            case MaterialLib::Solids::Phasefield::PhaseFieldModel::AT1:
            {
                element_surface_energy +=
                    gc * 0.375 *
                    ((1 - d_ip) / ls + (dNdx * d).dot((dNdx * d)) * ls) * w;

                break;
            }
            case MaterialLib::Solids::Phasefield::PhaseFieldModel::AT2:
            {
                element_surface_energy += 0.5 * gc *
                                          ((1 - d_ip) * (1 - d_ip) / ls +
                                           (dNdx * d).dot((dNdx * d)) * ls) *
                                          w;
                break;
            }
            case MaterialLib::Solids::Phasefield::PhaseFieldModel::COHESIVE:
            {
                element_surface_energy +=
                    gc / std::numbers::pi *
                    ((1 - d_ip * d_ip) / ls + (dNdx * d).dot((dNdx * d)) * ls) *
                    w;
                break;
            }
        }

        element_pressure_work += p_ip * (N_u_op(N) * u).dot(dNdx * d) * w;
    }

#ifdef USE_PETSC
    int const n_all_nodes = indices_of_processes[1].size();
    int const n_regular_nodes = std::count_if(
        begin(indices_of_processes[1]), end(indices_of_processes[1]),
        [](GlobalIndexType const& index) { return index >= 0; });
    if (n_all_nodes != n_regular_nodes)
    {
        element_elastic_energy *=
            static_cast<double>(n_regular_nodes) / n_all_nodes;
        element_surface_energy *=
            static_cast<double>(n_regular_nodes) / n_all_nodes;
        element_pressure_work *=
            static_cast<double>(n_regular_nodes) / n_all_nodes;
    }
#endif  // USE_PETSC
    elastic_energy += element_elastic_energy;
    surface_energy += element_surface_energy;
    pressure_work += element_pressure_work;
}

template <typename ShapeFunction, int DisplacementDim>
std::vector<double> const&
HMPhaseFieldLocalAssembler<ShapeFunction, DisplacementDim>::getIntPtWidth(
    const double /*t*/,
    std::vector<GlobalVector*> const& /*x*/,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
    std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(_ip_data,
                                                     &IpData::width_ip, cache);
}
}  // namespace HMPhaseField
}  // namespace ProcessLib
