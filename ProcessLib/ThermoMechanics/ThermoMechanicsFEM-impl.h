/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on July 2, 2019, 2:12 PM
 */

#pragma once

#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "ThermoMechanicsFEM.h"

namespace ProcessLib
{
namespace ThermoMechanics
{
template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
ThermoMechanicsLocalAssembler<ShapeFunction, IntegrationMethod,
                              DisplacementDim>::
    ThermoMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        ThermoMechanicsProcessData<DisplacementDim>& process_data)
    : _process_data(process_data),
      _integration_method(integration_order),
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
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
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
        ip_data.solid_density =
            _process_data.reference_solid_density(0, x_position)[0];
        ip_data.solid_density_prev = ip_data.solid_density;
        ip_data.N = shape_matrices[ip].N;
        ip_data.dNdx = shape_matrices[ip].dNdx;

        _secondary_data.N[ip] = shape_matrices[ip].N;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
std::size_t ThermoMechanicsLocalAssembler<
    ShapeFunction, IntegrationMethod,
    DisplacementDim>::setIPDataInitialConditions(std::string const& name,
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

    if (name == "sigma_ip")
    {
        return setSigma(values);
    }
    if (name == "epsilon_ip")
    {
        return setEpsilon(values);
    }
    if (name == "epsilon_m_ip")
    {
        return setEpsilonMechanical(values);
    }

    return 0;
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void ThermoMechanicsLocalAssembler<ShapeFunction, IntegrationMethod,
                                   DisplacementDim>::
    assembleWithJacobian(double const t, double const dt,
                         std::vector<double> const& local_x,
                         std::vector<double> const& local_xdot,
                         const double /*dxdot_dx*/, const double /*dx_dx*/,
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

    auto T_dot = Eigen::Map<typename ShapeMatricesType::template VectorType<
        temperature_size> const>(local_xdot.data() + temperature_index,
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

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

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

        double const dT = N.dot(T_dot) * dt;
        // calculate thermally induced strain
        // assume isotropic thermal expansion
        auto const alpha = _process_data.linear_thermal_expansion_coefficient(
            t, x_position)[0];
        double const linear_thermal_strain_increment = alpha * dT;

        //
        // displacement equation, displacement part
        //
        eps.noalias() = B * u;

        using Invariants = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value>;

        // assume isotropic thermal expansion
        const double T_ip = N.dot(T);  // T at integration point
        eps_m.noalias() =
            eps_m_prev + eps - eps_prev -
            linear_thermal_strain_increment * Invariants::identity2;

        auto&& solution = _ip_data[ip].solid_material.integrateStress(
            t, x_position, dt, eps_m_prev, eps_m, sigma_prev, *state, T_ip);

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

        // calculate real density
        // rho_s_{n+1} * (V_{n} + dV) = rho_s_n * V_n
        // dV = 3 * alpha * dT * V_0
        // rho_s_{n+1} = rho_s_n / (1 + 3 * alpha * dT )
        // see reference solid density description for details.
        auto& rho_s = _ip_data[ip].solid_density;
        rho_s = _ip_data[ip].solid_density_prev /
                (1 + 3 * linear_thermal_strain_increment);

        auto const& b = _process_data.specific_body_force;
        local_rhs.template block<displacement_size, 1>(displacement_index, 0)
            .noalias() -=
            (B.transpose() * sigma - N_u.transpose() * rho_s * b) * w;

        //
        // displacement equation, temperature part
        //
        KuT.noalias() +=
            B.transpose() * C * alpha * Invariants::identity2 * N * w;
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
            _process_data.thermal_conductivity(t, x_position)[0];
        KTT.noalias() += dNdx.transpose() * lambda * dNdx * w;

        auto const c = _process_data.specific_heat_capacity(t, x_position)[0];
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
        .noalias() -= KTT * T + DTT * T_dot;
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void ThermoMechanicsLocalAssembler<ShapeFunction, IntegrationMethod,
                                   DisplacementDim>::
    assembleWithJacobianForStaggeredScheme(
        const double t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_xdot, const double /*dxdot_dx*/,
        const double /*dx_dx*/, int const process_id,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    // For the equations with pressure
    if (process_id == _process_data.heat_conduction_process_id)
    {
        assembleWithJacobianForHeatConductionEquations(
            t, dt, local_x, local_xdot, local_b_data, local_Jac_data);
        return;
    }

    // For the equations with deformation
    assembleWithJacobianForDeformationEquations(t, dt, local_x, local_xdot,
                                                local_b_data, local_Jac_data);
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void ThermoMechanicsLocalAssembler<ShapeFunction, IntegrationMethod,
                                   DisplacementDim>::
    assembleWithJacobianForDeformationEquations(
        const double t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_xdot, std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data)
{
    auto const local_T =
        local_x.template segment<temperature_size>(temperature_index);

    auto const local_Tdot =
        local_xdot.template segment<temperature_size>(temperature_index);

    auto const local_u =
        local_x.template segment<displacement_size>(displacement_index);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesType::template MatrixType<displacement_size,
                                                        displacement_size>>(
        local_Jac_data, displacement_size, displacement_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesType::template VectorType<displacement_size>>(
        local_b_data, displacement_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

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
        double const dT_ip = N.dot(local_Tdot) * dt;
        // calculate thermally induced strain
        // assume isotropic thermal expansion
        auto const alpha = _process_data.linear_thermal_expansion_coefficient(
            t, x_position)[0];
        double const linear_thermal_strain_increment = alpha * dT_ip;

        //
        // displacement equation, displacement part
        //
        eps.noalias() = B * local_u;

        using Invariants = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value>;

        // assume isotropic thermal expansion
        eps_m.noalias() =
            eps_m_prev + eps - eps_prev -
            linear_thermal_strain_increment * Invariants::identity2;

        auto&& solution = _ip_data[ip].solid_material.integrateStress(
            t, x_position, dt, eps_m_prev, eps_m, sigma_prev, *state, T_ip);

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

        // calculate real density
        // rho_s_{n+1} * (V_{n} + dV) = rho_s_n * V_n
        // dV = 3 * alpha * dT * V_0
        // rho_s_{n+1} = rho_s_n / (1 + 3 * alpha * dT )
        // see reference solid density description for details.
        auto& rho_s = _ip_data[ip].solid_density;
        rho_s = _ip_data[ip].solid_density_prev /
                (1 + 3 * linear_thermal_strain_increment);

        auto const& b = _process_data.specific_body_force;
        local_rhs.noalias() -=
            (B.transpose() * sigma - N_u.transpose() * rho_s * b) * w;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void ThermoMechanicsLocalAssembler<ShapeFunction, IntegrationMethod,
                                   DisplacementDim>::
    assembleWithJacobianForHeatConductionEquations(
        const double t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_xdot, std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data)
{
    auto const local_T =
        local_x.template segment<temperature_size>(temperature_index);

    auto const local_dT =
        local_xdot.template segment<temperature_size>(temperature_index) * dt;

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

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        // calculate real density
        // rho_s_{n+1} * (V_{n} + dV) = rho_s_n * V_n
        // dV = 3 * alpha * dT * V_0
        // rho_s_{n+1} = rho_s_n / (1 + 3 * alpha * dT )
        // see reference solid density description for details.
        auto& rho_s = _ip_data[ip].solid_density;
        // calculate thermally induced strain
        // assume isotropic thermal expansion
        auto const alpha = _process_data.linear_thermal_expansion_coefficient(
            t, x_position)[0];

        double const dT_ip = N.dot(local_dT);
        double const linear_thermal_strain_increment = alpha * dT_ip;
        rho_s = _ip_data[ip].solid_density_prev /
                (1 + 3 * linear_thermal_strain_increment);
        auto const c_p = _process_data.specific_heat_capacity(t, x_position)[0];
        mass.noalias() += N.transpose() * rho_s * c_p * N * w;

        auto const lambda =
            _process_data.thermal_conductivity(t, x_position)[0];
        laplace.noalias() += dNdx.transpose() * lambda * dNdx * w;
    }
    local_Jac.noalias() += laplace + mass / dt;

    local_rhs.noalias() -= laplace * local_T + mass * local_dT / dt;
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
std::size_t
ThermoMechanicsLocalAssembler<ShapeFunction, IntegrationMethod,
                              DisplacementDim>::setSigma(double const* values)
{
    return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
        values, _ip_data, &IpData::sigma);
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
std::vector<double> ThermoMechanicsLocalAssembler<
    ShapeFunction, IntegrationMethod, DisplacementDim>::getSigma() const
{
    return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
        _ip_data, &IpData::sigma);
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
std::vector<double> const& ThermoMechanicsLocalAssembler<
    ShapeFunction, IntegrationMethod, DisplacementDim>::
    getIntPtSigma(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
        _ip_data, &IpData::sigma, cache);
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
std::size_t
ThermoMechanicsLocalAssembler<ShapeFunction, IntegrationMethod,
                              DisplacementDim>::setEpsilon(double const* values)
{
    return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
        values, _ip_data, &IpData::eps);
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
std::vector<double> ThermoMechanicsLocalAssembler<
    ShapeFunction, IntegrationMethod, DisplacementDim>::getEpsilon() const
{
    return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
        _ip_data, &IpData::eps);
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
std::vector<double> const& ThermoMechanicsLocalAssembler<
    ShapeFunction, IntegrationMethod, DisplacementDim>::
    getIntPtEpsilon(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
        _ip_data, &IpData::eps, cache);
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
std::size_t ThermoMechanicsLocalAssembler<
    ShapeFunction, IntegrationMethod,
    DisplacementDim>::setEpsilonMechanical(double const* values)
{
    auto const kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    auto epsilon_m_values =
        Eigen::Map<Eigen::Matrix<double, kelvin_vector_size, Eigen::Dynamic,
                                 Eigen::ColMajor> const>(
            values, kelvin_vector_size, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        _ip_data[ip].eps_m =
            MathLib::KelvinVector::symmetricTensorToKelvinVector(
                epsilon_m_values.col(ip));
    }

    return n_integration_points;
}
template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
std::vector<double>
ThermoMechanicsLocalAssembler<ShapeFunction, IntegrationMethod,
                              DisplacementDim>::getEpsilonMechanical() const
{
    auto const kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    std::vector<double> ip_epsilon_m_values;
    auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, Eigen::Dynamic, kelvin_vector_size, Eigen::RowMajor>>(
        ip_epsilon_m_values, n_integration_points, kelvin_vector_size);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        auto const& eps_m = _ip_data[ip].eps_m;
        cache_mat.row(ip) =
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(eps_m);
    }

    return ip_epsilon_m_values;
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
unsigned ThermoMechanicsLocalAssembler<
    ShapeFunction, IntegrationMethod,
    DisplacementDim>::getNumberOfIntegrationPoints() const
{
    return _integration_method.getNumberOfPoints();
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
typename MaterialLib::Solids::MechanicsBase<
    DisplacementDim>::MaterialStateVariables const&
ThermoMechanicsLocalAssembler<ShapeFunction, IntegrationMethod,
                              DisplacementDim>::
    getMaterialStateVariablesAt(unsigned integration_point) const
{
    return *_ip_data[integration_point].material_state_variables;
}
}  // namespace ThermoMechanics
}  // namespace ProcessLib
