/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * Created on July 2, 2019, 2:12 PM
 */

#pragma once

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/MPL/Utils/FormKelvinVector.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "ProcessLib/Utils/TransposeInPlace.h"

namespace ProcessLib
{
namespace ThermoMechanics
{
template <typename ShapeFunction, int DisplacementDim>
ThermoMechanicsLocalAssembler<ShapeFunction, DisplacementDim>::
    ThermoMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        ThermoMechanicsProcessData<DisplacementDim>& process_data)
    : _process_data(process_data),
      _integration_method(integration_method),
      _element(e),
      _is_axially_symmetric(is_axially_symmetric)
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);
    _secondary_data.N.resize(n_integration_points);

    auto const shape_matrices =
        NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                  DisplacementDim>(e, is_axially_symmetric,
                                                   _integration_method);

    auto& solid_material = MaterialLib::Solids::selectSolidConstitutiveRelation(
        _process_data.solid_materials, _process_data.material_ids, e.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        _ip_data.emplace_back(solid_material);
        auto& ip_data = _ip_data[ip];
        ip_data.integration_weight =
            _integration_method.getWeightedPoint(ip).getWeight() *
            shape_matrices[ip].integralMeasure * shape_matrices[ip].detJ;

        static const int kelvin_vector_size =
            MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
        // Initialize current time step values
        ip_data.sigma.setZero(kelvin_vector_size);
        ip_data.eps.setZero(kelvin_vector_size);

        // Previous time step values are not initialized and are set later.
        ip_data.sigma_prev.resize(kelvin_vector_size);
        ip_data.eps_prev.resize(kelvin_vector_size);

        ip_data.eps_m.setZero(kelvin_vector_size);
        ip_data.eps_m_prev.setZero(kelvin_vector_size);
        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(_element.getID());
        ip_data.N = shape_matrices[ip].N;
        ip_data.dNdx = shape_matrices[ip].dNdx;

        _secondary_data.N[ip] = shape_matrices[ip].N;
    }
}

template <typename ShapeFunction, int DisplacementDim>
std::size_t ThermoMechanicsLocalAssembler<ShapeFunction, DisplacementDim>::
    setIPDataInitialConditions(std::string_view const name,
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
        return setSigma(values);
    }
    if (name == "epsilon")
    {
        return setEpsilon(values);
    }
    if (name == "epsilon_m")
    {
        return setEpsilonMechanical(values);
    }

    return 0;
}

template <typename ShapeFunction, int DisplacementDim>
void ThermoMechanicsLocalAssembler<ShapeFunction, DisplacementDim>::
    assembleWithJacobian(double const t, double const dt,
                         std::vector<double> const& local_x,
                         std::vector<double> const& local_x_prev,
                         std::vector<double>& /*local_M_data*/,
                         std::vector<double>& /*local_K_data*/,
                         std::vector<double>& local_rhs_data,
                         std::vector<double>& local_Jac_data)
{
    auto const local_matrix_size = local_x.size();
    assert(local_matrix_size == temperature_size + displacement_size);

    auto T = Eigen::Map<typename ShapeMatricesType::template VectorType<
        temperature_size> const>(local_x.data() + temperature_index,
                                 temperature_size);

    auto u = Eigen::Map<typename ShapeMatricesType::template VectorType<
        displacement_size> const>(local_x.data() + displacement_index,
                                  displacement_size);
    bool const is_u_non_zero = u.norm() > 0.0;

    auto T_prev = Eigen::Map<typename ShapeMatricesType::template VectorType<
        temperature_size> const>(local_x_prev.data() + temperature_index,
                                 temperature_size);

    auto local_Jac = MathLib::createZeroedMatrix<JacobianMatrix>(
        local_Jac_data, local_matrix_size, local_matrix_size);

    auto local_rhs = MathLib::createZeroedVector<RhsVector>(local_rhs_data,
                                                            local_matrix_size);

    typename ShapeMatricesType::template MatrixType<displacement_size,
                                                    temperature_size>
        KuT;
    KuT.setZero(displacement_size, temperature_size);

    typename ShapeMatricesType::NodalMatrixType KTT;
    KTT.setZero(temperature_size, temperature_size);

    typename ShapeMatricesType::NodalMatrixType DTT;
    DTT.setZero(temperature_size, temperature_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;
    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    auto const& medium = _process_data.media_map.getMedium(_element.getID());
    auto const& solid_phase = medium->phase("Solid");

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                _element, N);
        auto const& B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunction::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx, N, x_coord, _is_axially_symmetric);

        auto& sigma = _ip_data[ip].sigma;
        auto const& sigma_prev = _ip_data[ip].sigma_prev;

        auto& eps = _ip_data[ip].eps;
        auto const& eps_prev = _ip_data[ip].eps_prev;

        auto& eps_m = _ip_data[ip].eps_m;
        auto const& eps_m_prev = _ip_data[ip].eps_m_prev;

        auto& state = _ip_data[ip].material_state_variables;

        const double T_ip = N.dot(T);  // T at integration point
        double const T_prev_ip = N.dot(T_prev);

        // Consider also anisotropic thermal expansion.
        auto const solid_linear_thermal_expansivity_vector =
            MPL::formKelvinVector<DisplacementDim>(
                solid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_expansivity)
                    .value(variables, x_position, t, dt));

        MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const
            dthermal_strain =
                solid_linear_thermal_expansivity_vector * (T_ip - T_prev_ip);

        //
        // displacement equation, displacement part
        //
        // For the restart computation, the displacement may not be
        // reloaded but the initial strains are always available. For such case,
        // the following computation is skipped.
        if (is_u_non_zero)
        {
            eps.noalias() = B * u;
        }

        eps_m.noalias() = eps_m_prev + eps - eps_prev - dthermal_strain;

        variables_prev.stress
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                sigma_prev);
        variables_prev.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m_prev);
        variables_prev.temperature = T_ip;
        variables.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m);
        variables.temperature = T_ip;

        auto&& solution = _ip_data[ip].solid_material.integrateStress(
            variables_prev, variables, t, x_position, dt, *state);

        if (!solution)
        {
            OGS_FATAL("Computation of local constitutive relation failed.");
        }

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
        std::tie(sigma, state, C) = std::move(*solution);

        local_Jac
            .template block<displacement_size, displacement_size>(
                displacement_index, displacement_index)
            .noalias() += B.transpose() * C * B * w;

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

        auto const rho_s =
            solid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);

        auto const& b = _process_data.specific_body_force;
        local_rhs.template block<displacement_size, 1>(displacement_index, 0)
            .noalias() -=
            (B.transpose() * sigma - N_u.transpose() * rho_s * b) * w;

        //
        // displacement equation, temperature part
        // The computation of KuT can be ignored.
        auto const alpha_T_tensor =
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                solid_linear_thermal_expansivity_vector);
        KuT.noalias() += B.transpose() * (C * alpha_T_tensor) * N * w;

        if (_ip_data[ip].solid_material.getConstitutiveModel() ==
            MaterialLib::Solids::ConstitutiveModel::CreepBGRa)
        {
            auto const s = Invariants::deviatoric_projection * sigma;
            double const norm_s = Invariants::FrobeniusNorm(s);
            const double creep_coefficient =
                _ip_data[ip].solid_material.getTemperatureRelatedCoefficient(
                    t, dt, x_position, T_ip, norm_s);
            KuT.noalias() += creep_coefficient * B.transpose() * s * N * w;
        }

        //
        // temperature equation, temperature part;
        //
        auto const lambda =
            solid_phase
                .property(
                    MaterialPropertyLib::PropertyType::thermal_conductivity)
                .value(variables, x_position, t, dt);

        GlobalDimMatrixType const thermal_conductivity =
            MaterialPropertyLib::formEigenTensor<DisplacementDim>(lambda);

        KTT.noalias() += dNdx.transpose() * thermal_conductivity * dNdx * w;

        auto const c =
            solid_phase
                .property(
                    MaterialPropertyLib::PropertyType::specific_heat_capacity)
                .template value<double>(variables, x_position, t, dt);
        DTT.noalias() += N.transpose() * rho_s * c * N * w;
    }

    // temperature equation, temperature part
    local_Jac
        .template block<temperature_size, temperature_size>(temperature_index,
                                                            temperature_index)
        .noalias() += KTT + DTT / dt;

    // displacement equation, temperature part
    local_Jac
        .template block<displacement_size, temperature_size>(displacement_index,
                                                             temperature_index)
        .noalias() -= KuT;

    local_rhs.template block<temperature_size, 1>(temperature_index, 0)
        .noalias() -= KTT * T + DTT * (T - T_prev) / dt;
}

template <typename ShapeFunction, int DisplacementDim>
void ThermoMechanicsLocalAssembler<ShapeFunction, DisplacementDim>::
    assembleWithJacobianForStaggeredScheme(
        const double t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_x_prev, int const process_id,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    // For the equations with pressure
    if (process_id == _process_data.heat_conduction_process_id)
    {
        assembleWithJacobianForHeatConductionEquations(
            t, dt, local_x, local_x_prev, local_b_data, local_Jac_data);
        return;
    }

    // For the equations with deformation
    assembleWithJacobianForDeformationEquations(t, dt, local_x, local_x_prev,
                                                local_b_data, local_Jac_data);
}

template <typename ShapeFunction, int DisplacementDim>
void ThermoMechanicsLocalAssembler<ShapeFunction, DisplacementDim>::
    assembleWithJacobianForDeformationEquations(
        const double t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_x_prev, std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data)
{
    auto const local_T =
        local_x.template segment<temperature_size>(temperature_index);

    auto const local_T_prev =
        local_x_prev.template segment<temperature_size>(temperature_index);

    auto const local_u =
        local_x.template segment<displacement_size>(displacement_index);
    bool const is_u_non_zero = local_u.norm() > 0.0;

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesType::template MatrixType<displacement_size,
                                                        displacement_size>>(
        local_Jac_data, displacement_size, displacement_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesType::template VectorType<displacement_size>>(
        local_b_data, displacement_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;
    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());
    auto const& medium = _process_data.media_map.getMedium(_element.getID());
    auto const& solid_phase = medium->phase("Solid");

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                _element, N);
        auto const& B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunction::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx, N, x_coord, _is_axially_symmetric);

        auto& sigma = _ip_data[ip].sigma;
        auto const& sigma_prev = _ip_data[ip].sigma_prev;

        auto& eps = _ip_data[ip].eps;
        auto const& eps_prev = _ip_data[ip].eps_prev;

        auto& eps_m = _ip_data[ip].eps_m;
        auto const& eps_m_prev = _ip_data[ip].eps_m_prev;

        auto& state = _ip_data[ip].material_state_variables;

        const double T_ip = N.dot(local_T);  // T at integration point
        variables.temperature = T_ip;
        double const dT_ip = T_ip - N.dot(local_T_prev);

        //
        // displacement equation, displacement part
        //
        // For the restart computation, the displacement may not be
        // reloaded but the initial strains are always available. For such case,
        // the following computation is skipped.
        if (is_u_non_zero)
        {
            eps.noalias() = B * local_u;
        }

        // Consider also anisotropic thermal expansion.
        auto const solid_linear_thermal_expansivity_vector =
            MPL::formKelvinVector<DisplacementDim>(
                solid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_expansivity)
                    .value(variables, x_position, t, dt));

        MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const
            dthermal_strain = solid_linear_thermal_expansivity_vector * dT_ip;

        eps_m.noalias() = eps_m_prev + eps - eps_prev - dthermal_strain;

        variables_prev.stress
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                sigma_prev);
        variables_prev.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m_prev);
        variables_prev.temperature = T_ip;
        variables.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m);

        auto&& solution = _ip_data[ip].solid_material.integrateStress(
            variables_prev, variables, t, x_position, dt, *state);

        if (!solution)
        {
            OGS_FATAL("Computation of local constitutive relation failed.");
        }

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
        std::tie(sigma, state, C) = std::move(*solution);

        local_Jac.noalias() += B.transpose() * C * B * w;

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

        auto const rho_s =
            solid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);

        auto const& b = _process_data.specific_body_force;
        local_rhs.noalias() -=
            (B.transpose() * sigma - N_u.transpose() * rho_s * b) * w;
    }
}

template <typename ShapeFunction, int DisplacementDim>
void ThermoMechanicsLocalAssembler<ShapeFunction, DisplacementDim>::
    assembleWithJacobianForHeatConductionEquations(
        const double t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_x_prev, std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data)
{
    auto const local_T =
        local_x.template segment<temperature_size>(temperature_index);

    auto const local_T_prev =
        local_x_prev.template segment<temperature_size>(temperature_index);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesType::template MatrixType<temperature_size,
                                                        temperature_size>>(
        local_Jac_data, temperature_size, temperature_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesType::template VectorType<temperature_size>>(
        local_b_data, temperature_size);

    typename ShapeMatricesType::NodalMatrixType mass;
    mass.setZero(temperature_size, temperature_size);

    typename ShapeMatricesType::NodalMatrixType laplace;
    laplace.setZero(temperature_size, temperature_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());
    auto const& medium = _process_data.media_map.getMedium(_element.getID());
    auto const& solid_phase = medium->phase("Solid");
    MPL::VariableArray variables;

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        const double T_ip = N.dot(local_T);  // T at integration point
        variables.temperature = T_ip;

        auto const rho_s =
            solid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);
        auto const c_p =
            solid_phase.property(MPL::PropertyType::specific_heat_capacity)
                .template value<double>(variables, x_position, t, dt);

        mass.noalias() += N.transpose() * rho_s * c_p * N * w;

        auto const lambda =
            solid_phase
                .property(
                    MaterialPropertyLib::PropertyType::thermal_conductivity)
                .value(variables, x_position, t, dt);

        GlobalDimMatrixType const thermal_conductivity =
            MaterialPropertyLib::formEigenTensor<DisplacementDim>(lambda);

        laplace.noalias() += dNdx.transpose() * thermal_conductivity * dNdx * w;
    }
    local_Jac.noalias() += laplace + mass / dt;

    local_rhs.noalias() -=
        laplace * local_T + mass * (local_T - local_T_prev) / dt;
}

template <typename ShapeFunction, int DisplacementDim>
std::size_t
ThermoMechanicsLocalAssembler<ShapeFunction, DisplacementDim>::setSigma(
    double const* values)
{
    return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
        values, _ip_data, &IpData::sigma);
}

template <typename ShapeFunction, int DisplacementDim>
std::vector<double>
ThermoMechanicsLocalAssembler<ShapeFunction, DisplacementDim>::getSigma() const
{
    constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);

    return transposeInPlace<kelvin_vector_size>(
        [this](std::vector<double>& values)
        { return getIntPtSigma(0, {}, {}, values); });
}

template <typename ShapeFunction, int DisplacementDim>
std::vector<double> const&
ThermoMechanicsLocalAssembler<ShapeFunction, DisplacementDim>::getIntPtSigma(
    const double /*t*/,
    std::vector<GlobalVector*> const& /*x*/,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
    std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
        _ip_data, &IpData::sigma, cache);
}

template <typename ShapeFunction, int DisplacementDim>
std::size_t
ThermoMechanicsLocalAssembler<ShapeFunction, DisplacementDim>::setEpsilon(
    double const* values)
{
    return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
        values, _ip_data, &IpData::eps);
}

template <typename ShapeFunction, int DisplacementDim>
std::vector<double> ThermoMechanicsLocalAssembler<
    ShapeFunction, DisplacementDim>::getEpsilon() const
{
    constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);

    return transposeInPlace<kelvin_vector_size>(
        [this](std::vector<double>& values)
        { return getIntPtEpsilon(0, {}, {}, values); });
}

template <typename ShapeFunction, int DisplacementDim>
std::vector<double> const&
ThermoMechanicsLocalAssembler<ShapeFunction, DisplacementDim>::getIntPtEpsilon(
    const double /*t*/,
    std::vector<GlobalVector*> const& /*x*/,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
    std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
        _ip_data, &IpData::eps, cache);
}

template <typename ShapeFunction, int DisplacementDim>
std::vector<double> const&
ThermoMechanicsLocalAssembler<ShapeFunction, DisplacementDim>::
    getIntPtEpsilonMechanical(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
        _ip_data, &IpData::eps_m, cache);
}

template <typename ShapeFunction, int DisplacementDim>
std::size_t ThermoMechanicsLocalAssembler<
    ShapeFunction, DisplacementDim>::setEpsilonMechanical(double const* values)
{
    return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
        values, _ip_data, &IpData::eps_m);
}
template <typename ShapeFunction, int DisplacementDim>
std::vector<double> ThermoMechanicsLocalAssembler<
    ShapeFunction, DisplacementDim>::getEpsilonMechanical() const
{
    constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);

    return transposeInPlace<kelvin_vector_size>(
        [this](std::vector<double>& values)
        { return getIntPtEpsilonMechanical(0, {}, {}, values); });
}

template <typename ShapeFunction, int DisplacementDim>
unsigned ThermoMechanicsLocalAssembler<
    ShapeFunction, DisplacementDim>::getNumberOfIntegrationPoints() const
{
    return _integration_method.getNumberOfPoints();
}

template <typename ShapeFunction, int DisplacementDim>
int ThermoMechanicsLocalAssembler<ShapeFunction,
                                  DisplacementDim>::getMaterialID() const
{
    return _process_data.material_ids == nullptr
               ? 0
               : (*_process_data.material_ids)[_element.getID()];
}

template <typename ShapeFunction, int DisplacementDim>
typename MaterialLib::Solids::MechanicsBase<
    DisplacementDim>::MaterialStateVariables const&
ThermoMechanicsLocalAssembler<ShapeFunction, DisplacementDim>::
    getMaterialStateVariablesAt(unsigned integration_point) const
{
    return *_ip_data[integration_point].material_state_variables;
}
}  // namespace ThermoMechanics
}  // namespace ProcessLib
