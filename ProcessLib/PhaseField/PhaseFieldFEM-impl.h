// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <numbers>

#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
namespace PhaseField
{
template <typename ShapeFunction, int DisplacementDim>
void PhaseFieldLocalAssembler<ShapeFunction, DisplacementDim>::
    assembleWithJacobianForStaggeredScheme(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& /*local_x_prev*/, int const process_id,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    // For the equations with phase field.
    if (process_id == phase_process_id)
    {
        assembleWithJacobianPhaseFieldEquations(t, dt, local_x, local_b_data,
                                                local_Jac_data);
        return;
    }

    // For the equations with deformation
    assembleWithJacobianForDeformationEquations(t, dt, local_x, local_b_data,
                                                local_Jac_data);
}

template <typename ShapeFunction, int DisplacementDim>
void PhaseFieldLocalAssembler<ShapeFunction, DisplacementDim>::
    assembleWithJacobianForDeformationEquations(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    using DeformationMatrix =
        typename ShapeMatricesType::template MatrixType<displacement_size,
                                                        displacement_size>;
    auto const d = local_x.template segment<phasefield_size>(phasefield_index);
    auto const u =
        local_x.template segment<displacement_size>(displacement_index);

    auto local_Jac = MathLib::createZeroedMatrix<DeformationMatrix>(
        local_Jac_data, displacement_size, displacement_size);

    auto local_rhs = MathLib::createZeroedVector<DeformationVector>(
        local_b_data, displacement_size);

    auto local_pressure = 0.0;
    if (_process_data.pressurized_crack)
    {
        local_pressure = _process_data.unity_pressure;
    }

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

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
        double const d_ip = N.dot(d);
        double const degradation =
            _process_data.degradation_derivative->degradation(d_ip, k, ls);
        _ip_data[ip].updateConstitutiveRelation(
            t, x_position, dt, u, degradation,
            _process_data.energy_split_model);

        auto& sigma = _ip_data[ip].sigma;
        auto const& D = _ip_data[ip].D;

        typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                        displacement_size>
            N_u = ShapeMatricesType::template MatrixType<
                DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                          displacement_size);

        for (int i = 0; i < DisplacementDim; ++i)
        {
            N_u.template block<1, displacement_size / DisplacementDim>(
                   i, i * displacement_size / DisplacementDim)
                .noalias() = N;
        }

        auto const rho_sr = _process_data.solid_density(t, x_position)[0];
        auto const& b = _process_data.specific_body_force;

        local_rhs.noalias() -=
            (B.transpose() * sigma - N_u.transpose() * rho_sr * b -
             local_pressure * N_u.transpose() * dNdx * d) *
            w;
        local_Jac.noalias() += B.transpose() * D * B * w;
    }
}

template <typename ShapeFunction, int DisplacementDim>
void PhaseFieldLocalAssembler<ShapeFunction, DisplacementDim>::
    assembleWithJacobianPhaseFieldEquations(double const t, double const dt,
                                            Eigen::VectorXd const& local_x,
                                            std::vector<double>& local_b_data,
                                            std::vector<double>& local_Jac_data)
{
    auto const d = local_x.template segment<phasefield_size>(phasefield_index);
    auto const u =
        local_x.template segment<displacement_size>(displacement_index);

    auto local_Jac = MathLib::createZeroedMatrix<PhaseFieldMatrix>(
        local_Jac_data, phasefield_size, phasefield_size);
    auto local_rhs = MathLib::createZeroedVector<PhaseFieldVector>(
        local_b_data, phasefield_size);

    auto local_pressure = 0.0;
    if (_process_data.static_pressurized_crack)
    {
        local_pressure = _process_data.unity_pressure;
    }
    else if (_process_data.propagating_pressurized_crack)
    {
        local_pressure = _process_data.pressure;
    }

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        double const d_ip = N.dot(d);

        ParameterLib::SpatialPosition const x_position{
            std::nullopt, this->_element.getID(),
            MathLib::Point3d(NumLib::interpolateCoordinates<ShapeFunction,
                                                            ShapeMatricesType>(
                this->_element, N))};

        double const k = _process_data.residual_stiffness(t, x_position)[0];
        double const ls = _process_data.crack_length_scale(t, x_position)[0];
        double const gc = _process_data.crack_resistance(t, x_position)[0];
        double const degradation =
            _process_data.degradation_derivative->degradation(d_ip, k, ls);
        double const degradation_df1 =
            _process_data.degradation_derivative->degradationDf1(d_ip, k, ls);
        double const degradation_df2 =
            _process_data.degradation_derivative->degradationDf2(d_ip, k, ls);

        auto const x_coord =
            x_position.getCoordinates().value()[0];  // r for axisymmetry
        auto const& B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunction::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx, N, x_coord, _is_axially_symmetric);

        auto& eps = _ip_data[ip].eps;
        eps.noalias() = B * u;
        _ip_data[ip].updateConstitutiveRelation(
            t, x_position, dt, u, degradation,
            _process_data.energy_split_model);

        auto const& strain_energy_tensile = _ip_data[ip].strain_energy_tensile;

        auto& ip_data = _ip_data[ip];
        ip_data.strain_energy_tensile = strain_energy_tensile;

        typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                        displacement_size>
            N_u = ShapeMatricesType::template MatrixType<
                DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                          displacement_size);

        for (int i = 0; i < DisplacementDim; ++i)
        {
            N_u.template block<1, displacement_size / DisplacementDim>(
                   i, i * displacement_size / DisplacementDim)
                .noalias() = N;
        }

        local_Jac.noalias() +=
            (N.transpose() * N * degradation_df2 * strain_energy_tensile) * w;

        local_rhs.noalias() -=
            (N.transpose() * degradation_df1 * strain_energy_tensile -
             local_pressure * dNdx.transpose() * N_u * u) *
            w;

        calculateCrackLocalJacobianAndResidual<
            decltype(dNdx), decltype(N), decltype(w), decltype(d),
            decltype(local_Jac), decltype(local_rhs)>(
            dNdx, N, w, d, local_Jac, local_rhs, gc, ls,
            _process_data.phasefield_model);
    }
}

template <typename ShapeFunction, int DisplacementDim>
void PhaseFieldLocalAssembler<ShapeFunction, DisplacementDim>::
    computeCrackIntegral(
        std::size_t mesh_item_id,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
        std::vector<GlobalVector*> const& x, double const /*t*/,
        double& crack_volume)
{
    std::vector<std::vector<GlobalIndexType>> indices_of_processes;
    indices_of_processes.reserve(dof_tables.size());
    std::transform(dof_tables.begin(), dof_tables.end(),
                   std::back_inserter(indices_of_processes),
                   [&](auto const dof_table)
                   { return NumLib::getIndices(mesh_item_id, *dof_table); });

    auto local_coupled_xs = getCoupledLocalSolutions(x, indices_of_processes);
    assert(local_coupled_xs.size() == displacement_size + phasefield_size);

    auto const d = Eigen::Map<PhaseFieldVector const>(
        &local_coupled_xs[phasefield_index], phasefield_size);
    auto const u = Eigen::Map<DeformationVector const>(
        &local_coupled_xs[displacement_index], displacement_size);

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                        displacement_size>
            N_u = ShapeMatricesType::template MatrixType<
                DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                          displacement_size);

        for (int i = 0; i < DisplacementDim; ++i)
        {
            N_u.template block<1, displacement_size / DisplacementDim>(
                   i, i * displacement_size / DisplacementDim)
                .noalias() = N;
        }

        crack_volume += (N_u * u).dot(dNdx * d) * w;
    }
}

template <typename ShapeFunction, int DisplacementDim>
void PhaseFieldLocalAssembler<ShapeFunction, DisplacementDim>::computeEnergy(
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
    assert(local_coupled_xs.size() == displacement_size + phasefield_size);

    auto const d = Eigen::Map<PhaseFieldVector const>(
        &local_coupled_xs[phasefield_index], phasefield_size);
    auto const u = Eigen::Map<DeformationVector const>(
        &local_coupled_xs[displacement_index], displacement_size);

    double element_elastic_energy = 0.0;
    double element_surface_energy = 0.0;
    double element_pressure_work = 0.0;

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;
        auto pressure_ip = _process_data.pressure;

        ParameterLib::SpatialPosition const x_position{
            std::nullopt, this->_element.getID(),
            MathLib::Point3d(NumLib::interpolateCoordinates<ShapeFunction,
                                                            ShapeMatricesType>(
                this->_element, N))};

        typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                        displacement_size>
            N_u = ShapeMatricesType::template MatrixType<
                DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                          displacement_size);

        for (int i = 0; i < DisplacementDim; ++i)
        {
            N_u.template block<1, displacement_size / DisplacementDim>(
                   i, i * displacement_size / DisplacementDim)
                .noalias() = N;
        }

        element_elastic_energy += _ip_data[ip].elastic_energy * w;

        double const gc = _process_data.crack_resistance(t, x_position)[0];
        double const ls = _process_data.crack_length_scale(t, x_position)[0];
        double const d_ip = N.dot(d);
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

        if (_process_data.pressurized_crack)
        {
            element_pressure_work += pressure_ip * (N_u * u).dot(dNdx * d) * w;
        }
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

template <typename ShapeFunctionDisplacement, int DisplacementDim>
std::size_t PhaseFieldLocalAssembler<
    ShapeFunctionDisplacement,
    DisplacementDim>::setIPDataInitialConditions(std::string_view const name,
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

    if (name == "sigma")
    {
        if (_process_data.initial_stress.value != nullptr)
        {
            OGS_FATAL(
                "Setting initial conditions for stress from integration "
                "point data and from a parameter '{:s}' is not possible "
                "simultaneously.",
                _process_data.initial_stress.value->name);
        }

        return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
            values, _ip_data, &IpData::sigma);
    }

    return 0;
}

template <typename ShapeFunctionDisplacement, int DisplacementDim>
std::vector<double> PhaseFieldLocalAssembler<ShapeFunctionDisplacement,
                                             DisplacementDim>::getSigma() const
{
    return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
        _ip_data, &IpData::sigma);
}

template <typename ShapeFunctionDisplacement, int DisplacementDim>
std::vector<double> PhaseFieldLocalAssembler<
    ShapeFunctionDisplacement, DisplacementDim>::getEpsilon() const
{
    auto const kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
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

}  // namespace PhaseField
}  // namespace ProcessLib
