/**
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on November 29, 2017, 2:03 PM
 */

#pragma once

#include <cassert>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/MPL/Utils/FormKelvinVectorFromThermalExpansivity.h"
#include "MaterialLib/MPL/Utils/GetLiquidThermalExpansivity.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "ProcessLib/Utils/TransposeInPlace.h"

namespace ProcessLib
{
namespace ThermoRichardsMechanics
{
template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement, ShapeFunction,
                                      IntegrationMethod, DisplacementDim>::
    ThermoRichardsMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        ThermoRichardsMechanicsProcessData<DisplacementDim>& process_data)
    : process_data_(process_data),
      integration_method_(integration_order),
      element_(e),
      is_axially_symmetric_(is_axially_symmetric)
{
    unsigned const n_integration_points =
        integration_method_.getNumberOfPoints();

    ip_data_.reserve(n_integration_points);
    secondary_data_.N_u.resize(n_integration_points);

    auto const shape_matrices_u =
        NumLib::initShapeMatrices<ShapeFunctionDisplacement,
                                  ShapeMatricesTypeDisplacement,
                                  DisplacementDim>(e, is_axially_symmetric,
                                                   integration_method_);

    auto const shape_matrices =
        NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                  DisplacementDim>(e, is_axially_symmetric,
                                                   integration_method_);

    auto const& solid_material =
        MaterialLib::Solids::selectSolidConstitutiveRelation(
            process_data_.solid_materials, process_data_.material_ids,
            e.getID());

    auto const& medium = process_data_.media_map->getMedium(element_.getID());

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(element_.getID());
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        ip_data_.emplace_back(solid_material);
        auto& ip_data = ip_data_[ip];
        auto const& sm_u = shape_matrices_u[ip];
        ip_data_[ip].integration_weight =
            integration_method_.getWeightedPoint(ip).getWeight() *
            sm_u.integralMeasure * sm_u.detJ;

        ip_data.N_u_op = ShapeMatricesTypeDisplacement::template MatrixType<
            DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                      displacement_size);
        for (int i = 0; i < DisplacementDim; ++i)
        {
            ip_data.N_u_op
                .template block<1, displacement_size / DisplacementDim>(
                    i, i * displacement_size / DisplacementDim)
                .noalias() = sm_u.N;
        }

        ip_data.N_u = sm_u.N;
        ip_data.dNdx_u = sm_u.dNdx;

        // ip_data.N_p and ip_data.dNdx_p are used for both p and T variables
        ip_data.N_p = shape_matrices[ip].N;
        ip_data.dNdx_p = shape_matrices[ip].dNdx;

        // Initial porosity. Could be read from integration point data or mesh.
        ip_data.porosity =
            medium->property(MPL::porosity)
                .template initialValue<double>(
                    x_position,
                    std::numeric_limits<
                        double>::quiet_NaN() /* t independent */);

        ip_data.transport_porosity = ip_data.porosity;
        if (medium->hasProperty(MPL::PropertyType::transport_porosity))
        {
            ip_data.transport_porosity =
                medium->property(MPL::transport_porosity)
                    .template initialValue<double>(
                        x_position,
                        std::numeric_limits<
                            double>::quiet_NaN() /* t independent */);
        }

        secondary_data_.N_u[ip] = shape_matrices_u[ip].N;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
std::size_t ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunction, IntegrationMethod,
    DisplacementDim>::setIPDataInitialConditions(std::string const& name,
                                                 double const* values,
                                                 int const integration_order)
{
    if (integration_order !=
        static_cast<int>(integration_method_.getIntegrationOrder()))
    {
        OGS_FATAL(
            "Setting integration point initial conditions; The integration "
            "order of the local assembler for element {:d} is different "
            "from the integration order in the initial condition.",
            element_.getID());
    }

    if (name == "sigma_ip")
    {
        if (process_data_.initial_stress != nullptr)
        {
            OGS_FATAL(
                "Setting initial conditions for stress from integration "
                "point data and from a parameter '{:s}' is not possible "
                "simultaneously.",
                process_data_.initial_stress->name);
        }
        return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
            values, ip_data_, &IpData::sigma_eff);
    }

    if (name == "saturation_ip")
    {
        return ProcessLib::setIntegrationPointScalarData(values, ip_data_,
                                                         &IpData::saturation);
    }
    if (name == "porosity_ip")
    {
        return ProcessLib::setIntegrationPointScalarData(values, ip_data_,
                                                         &IpData::porosity);
    }
    if (name == "transport_porosity_ip")
    {
        return ProcessLib::setIntegrationPointScalarData(
            values, ip_data_, &IpData::transport_porosity);
    }
    if (name == "swelling_stress_ip")
    {
        return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
            values, ip_data_, &IpData::sigma_sw);
    }
    if (name == "epsilon_ip")
    {
        return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
            values, ip_data_, &IpData::eps);
    }
    return 0;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
void ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                           ShapeFunction, IntegrationMethod,
                                           DisplacementDim>::
    setInitialConditionsConcrete(std::vector<double> const& local_x,
                                 double const t,
                                 bool const /*use_monolithic_scheme*/,
                                 int const /*process_id*/)
{
    assert(local_x.size() ==
           temperature_size + pressure_size + displacement_size);

    auto const p_L = Eigen::Map<
        typename ShapeMatricesType::template VectorType<pressure_size> const>(
        local_x.data() + pressure_index, pressure_size);

    auto const T = Eigen::Map<typename ShapeMatricesType::template VectorType<
        temperature_size> const>(local_x.data() + temperature_index,
                                 temperature_size);

    constexpr double dt = std::numeric_limits<double>::quiet_NaN();
    auto const& medium = process_data_.media_map->getMedium(element_.getID());
    MPL::VariableArray variables;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(element_.getID());

    auto const& solid_phase = medium->phase("Solid");

    unsigned const n_integration_points =
        integration_method_.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        // N is used for both T and p variables.
        auto const& N = ip_data_[ip].N_p;

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N, p_cap_ip);

        variables[static_cast<int>(MPL::Variable::capillary_pressure)] =
            p_cap_ip;
        variables[static_cast<int>(MPL::Variable::phase_pressure)] = -p_cap_ip;

        double T_ip;
        NumLib::shapeFunctionInterpolate(T, N, T_ip);
        variables[static_cast<int>(MPL::Variable::temperature)] = T_ip;

        ip_data_[ip].saturation_prev =
            medium->property(MPL::PropertyType::saturation)
                .template value<double>(variables, x_position, t, dt);

        // Set eps_m_prev from potentially non-zero eps and sigma_sw from
        // restart.
        auto const C_el = ip_data_[ip].computeElasticTangentStiffness(
            t, x_position, dt, T_ip, T_ip);
        auto& eps = ip_data_[ip].eps;
        auto& sigma_sw = ip_data_[ip].sigma_sw;
        ip_data_[ip].eps_m_prev.noalias() =
            solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate)
                ? eps + C_el.inverse() * sigma_sw
                : eps;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
void ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                           ShapeFunction, IntegrationMethod,
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
    auto const local_matrix_dim =
        displacement_size + pressure_size + temperature_size;
    assert(local_x.size() == local_matrix_dim);

    auto const T = Eigen::Map<typename ShapeMatricesType::template VectorType<
        temperature_size> const>(local_x.data() + temperature_index,
                                 temperature_size);
    auto const p_L = Eigen::Map<
        typename ShapeMatricesType::template VectorType<pressure_size> const>(
        local_x.data() + pressure_index, pressure_size);

    auto const u =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x.data() + displacement_index,
                                      displacement_size);

    auto const T_dot =
        Eigen::Map<typename ShapeMatricesType::template VectorType<
            temperature_size> const>(local_xdot.data() + temperature_index,
                                     temperature_size);
    auto const p_L_dot = Eigen::Map<
        typename ShapeMatricesType::template VectorType<pressure_size> const>(
        local_xdot.data() + pressure_index, pressure_size);
    auto const u_dot =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_xdot.data() + displacement_index,
                                      displacement_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            local_matrix_dim, local_matrix_dim>>(
        local_Jac_data, local_matrix_dim, local_matrix_dim);

    auto local_rhs =
        MathLib::createZeroedVector<typename ShapeMatricesTypeDisplacement::
                                        template VectorType<local_matrix_dim>>(
            local_rhs_data, local_matrix_dim);

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    typename ShapeMatricesType::NodalMatrixType M_TT =
        ShapeMatricesType::NodalMatrixType::Zero(temperature_size,
                                                 temperature_size);
    typename ShapeMatricesType::NodalMatrixType K_TT =
        ShapeMatricesType::NodalMatrixType::Zero(temperature_size,
                                                 temperature_size);
    typename ShapeMatricesType::NodalMatrixType K_Tp =
        ShapeMatricesType::NodalMatrixType::Zero(temperature_size,
                                                 pressure_size);
    typename ShapeMatricesType::NodalMatrixType M_pT =
        ShapeMatricesType::NodalMatrixType::Zero(pressure_size,
                                                 temperature_size);
    typename ShapeMatricesType::NodalMatrixType laplace_p =
        ShapeMatricesType::NodalMatrixType::Zero(pressure_size, pressure_size);

    typename ShapeMatricesType::NodalMatrixType storage_p_a_p =
        ShapeMatricesType::NodalMatrixType::Zero(pressure_size, pressure_size);

    typename ShapeMatricesType::NodalMatrixType storage_p_a_S_Jpp =
        ShapeMatricesType::NodalMatrixType::Zero(pressure_size, pressure_size);

    typename ShapeMatricesType::NodalMatrixType storage_p_a_S =
        ShapeMatricesType::NodalMatrixType::Zero(pressure_size, pressure_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        pressure_size, displacement_size>
        Kpu = ShapeMatricesTypeDisplacement::template MatrixType<
            pressure_size, displacement_size>::Zero(pressure_size,
                                                    displacement_size);

    auto const& medium = process_data_.media_map->getMedium(element_.getID());
    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& solid_phase = medium->phase("Solid");
    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(element_.getID());

    unsigned const n_integration_points =
        integration_method_.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto& ip_data = ip_data_[ip];
        auto const& w = ip_data.integration_weight;

        auto const& N_u_op = ip_data.N_u_op;

        auto const& N_u = ip_data.N_u;
        auto const& dNdx_u = ip_data.dNdx_u;

        // N and dNdx are used for both p and T variables
        auto const& N = ip_data.N_p;
        auto const& dNdx = ip_data.dNdx_p;

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                element_, N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, is_axially_symmetric_);

        double T_ip;
        NumLib::shapeFunctionInterpolate(T, N, T_ip);
        double T_dot_ip;
        NumLib::shapeFunctionInterpolate(T_dot, N, T_dot_ip);
        double const dT = T_dot_ip * dt;

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N, p_cap_ip);

        double p_cap_dot_ip;
        NumLib::shapeFunctionInterpolate(-p_L_dot, N, p_cap_dot_ip);

        variables[static_cast<int>(MPL::Variable::capillary_pressure)] =
            p_cap_ip;
        variables[static_cast<int>(MPL::Variable::phase_pressure)] = -p_cap_ip;
        variables[static_cast<int>(MPL::Variable::temperature)] = T_ip;

        auto& eps = ip_data.eps;
        auto& eps_m = ip_data.eps_m;
        eps.noalias() = B * u;
        auto const& sigma_eff = ip_data.sigma_eff;
        auto& S_L = ip_data.saturation;
        auto const S_L_prev = ip_data.saturation_prev;
        auto const alpha =
            medium->property(MPL::PropertyType::biot_coefficient)
                .template value<double>(variables, x_position, t, dt);

        double const T_ip_prev = T_ip - dT;
        auto const C_el = ip_data.computeElasticTangentStiffness(
            t, x_position, dt, T_ip_prev, T_ip);

        auto const beta_SR =
            (1 - alpha) /
            ip_data.solid_material.getBulkModulus(t, x_position, &C_el);
        variables[static_cast<int>(MPL::Variable::grain_compressibility)] =
            beta_SR;

        auto const rho_LR =
            liquid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);
        auto const& b = process_data_.specific_body_force;

        S_L = medium->property(MPL::PropertyType::saturation)
                  .template value<double>(variables, x_position, t, dt);
        variables[static_cast<int>(MPL::Variable::liquid_saturation)] = S_L;
        variables_prev[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L_prev;

        // tangent derivative for Jacobian
        double const dS_L_dp_cap =
            medium->property(MPL::PropertyType::saturation)
                .template dValue<double>(variables,
                                         MPL::Variable::capillary_pressure,
                                         x_position, t, dt);
        // secant derivative from time discretization for storage
        // use tangent, if secant is not available
        double const DeltaS_L_Deltap_cap =
            (p_cap_dot_ip == 0) ? dS_L_dp_cap
                                : (S_L - S_L_prev) / (dt * p_cap_dot_ip);

        auto const chi = [medium, x_position, t, dt](double const S_L)
        {
            MPL::VariableArray vs;
            vs[static_cast<int>(MPL::Variable::liquid_saturation)] = S_L;
            return medium->property(MPL::PropertyType::bishops_effective_stress)
                .template value<double>(vs, x_position, t, dt);
        };
        double const chi_S_L = chi(S_L);
        double const chi_S_L_prev = chi(S_L_prev);

        variables[static_cast<int>(MPL::Variable::effective_pore_pressure)] =
            -chi_S_L * p_cap_ip;
        variables_prev[static_cast<int>(
            MPL::Variable::effective_pore_pressure)] =
            -chi_S_L_prev * (p_cap_ip - p_cap_dot_ip * dt);

        // Set volumetric strain rate for the general case without swelling.
        variables[static_cast<int>(MPL::Variable::volumetric_strain)]
            .emplace<double>(Invariants::trace(eps));
        variables_prev[static_cast<int>(MPL::Variable::volumetric_strain)]
            .emplace<double>(Invariants::trace(B * (u - u_dot * dt)));

        auto& phi = ip_data.porosity;
        {  // Porosity update

            variables_prev[static_cast<int>(MPL::Variable::porosity)] =
                ip_data.porosity_prev;
            phi = medium->property(MPL::PropertyType::porosity)
                      .template value<double>(variables, variables_prev,
                                              x_position, t, dt);
            variables[static_cast<int>(MPL::Variable::porosity)] = phi;
        }

        if (alpha < phi)
        {
            OGS_FATAL(
                "ThermoRichardsMechanics: Biot-coefficient {} is smaller than "
                "porosity {} in element/integration point {}/{}.",
                alpha, phi, element_.getID(), ip);
        }

        // Swelling and possibly volumetric strain rate update.
        if (solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate))
        {
            auto& sigma_sw = ip_data.sigma_sw;
            auto const& sigma_sw_prev = ip_data.sigma_sw_prev;

            // If there is swelling, compute it. Update volumetric strain rate,
            // s.t. it corresponds to the mechanical part only.
            sigma_sw = sigma_sw_prev;

            using DimMatrix = Eigen::Matrix<double, 3, 3>;
            auto const sigma_sw_dot =
                MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
                    solid_phase
                        .property(MPL::PropertyType::swelling_stress_rate)
                        .template value<DimMatrix>(variables, variables_prev,
                                                   x_position, t, dt));
            sigma_sw += sigma_sw_dot * dt;

            // !!! Misusing volumetric strain for mechanical volumetric
            // strain just to update the transport porosity !!!
            std::get<double>(variables[static_cast<int>(
                MPL::Variable::volumetric_strain)]) +=
                identity2.transpose() * C_el.inverse() * sigma_sw;
            std::get<double>(variables_prev[static_cast<int>(
                MPL::Variable::volumetric_strain)]) +=
                identity2.transpose() * C_el.inverse() * sigma_sw_prev;
        }

        if (solid_phase.hasProperty(MPL::PropertyType::transport_porosity))
        {
            variables_prev[static_cast<int>(
                MPL::Variable::transport_porosity)] =
                ip_data.transport_porosity_prev;

            ip_data.transport_porosity =
                solid_phase.property(MPL::PropertyType::transport_porosity)
                    .template value<double>(variables, variables_prev,
                                            x_position, t, dt);
            variables[static_cast<int>(MPL::Variable::transport_porosity)] =
                ip_data.transport_porosity;
        }
        else
        {
            variables[static_cast<int>(MPL::Variable::transport_porosity)] =
                phi;
        }

        //
        // displacement equation, displacement part
        //

        // Consider also anisotropic thermal expansion.
        MathLib::KelvinVector::KelvinVectorType<
            DisplacementDim> const solid_linear_thermal_expansivity_vector =
            MPL::formKelvinVectorFromThermalExpansivity<DisplacementDim>(
                solid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_expansivity)
                    .value(variables, x_position, t, dt));

        MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const
            dthermal_strain = solid_linear_thermal_expansivity_vector * dT;

        auto& eps_prev = ip_data.eps_prev;
        auto& eps_m_prev = ip_data.eps_m_prev;
        eps_m.noalias() = eps_m_prev + eps - eps_prev - dthermal_strain;

        if (solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate))
        {
            eps_m.noalias() +=
                C_el.inverse() * (ip_data.sigma_sw - ip_data.sigma_sw_prev);
        }

        variables[static_cast<int>(
                      MaterialPropertyLib::Variable::mechanical_strain)]
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m);

        auto C = ip_data.updateConstitutiveRelation(variables, t, x_position,
                                                    dt, T_ip_prev);

        local_Jac
            .template block<displacement_size, displacement_size>(
                displacement_index, displacement_index)
            .noalias() += B.transpose() * C * B * w;

        double const p_FR = -chi_S_L * p_cap_ip;
        // p_SR
        variables[static_cast<int>(MPL::Variable::solid_grain_pressure)] =
            p_FR - Invariants::trace(sigma_eff) / (3 * (1 - phi));
        auto const rho_SR =
            solid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);

        double const rho = rho_SR * (1 - phi) + S_L * phi * rho_LR;

        auto const sigma_total =
            (sigma_eff + alpha * chi_S_L * identity2 * p_cap_ip).eval();

        local_rhs.template segment<displacement_size>(displacement_index)
            .noalias() -=
            (B.transpose() * sigma_total - N_u_op.transpose() * rho * b) * w;

        //
        // displacement equation, pressure part
        //
        auto const dchi_dS_L =
            medium->property(MPL::PropertyType::bishops_effective_stress)
                .template dValue<double>(variables,
                                         MPL::Variable::liquid_saturation,
                                         x_position, t, dt);
        local_Jac
            .template block<displacement_size, pressure_size>(
                displacement_index, pressure_index)
            .noalias() -= B.transpose() * alpha *
                          (chi_S_L + dchi_dS_L * p_cap_ip * dS_L_dp_cap) *
                          identity2 * N * w;

        local_Jac
            .template block<displacement_size, pressure_size>(
                displacement_index, pressure_index)
            .noalias() +=
            N_u_op.transpose() * phi * rho_LR * dS_L_dp_cap * b * N * w;

        if (solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate))
        {
            using DimMatrix = Eigen::Matrix<double, 3, 3>;
            auto const dsigma_sw_dS_L =
                MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
                    solid_phase
                        .property(MPL::PropertyType::swelling_stress_rate)
                        .template dValue<DimMatrix>(
                            variables, variables_prev,
                            MPL::Variable::liquid_saturation, x_position, t,
                            dt));
            local_Jac
                .template block<displacement_size, pressure_size>(
                    displacement_index, pressure_index)
                .noalias() +=
                B.transpose() * dsigma_sw_dS_L * dS_L_dp_cap * N * w;
        }

        //
        // pressure equation, displacement part.
        //
        Kpu.noalias() += N.transpose() * S_L * rho_LR * alpha *
                         identity2.transpose() * B * w;

        //
        // pressure equation, pressure part.
        //
        double const k_rel =
            medium->property(MPL::PropertyType::relative_permeability)
                .template value<double>(variables, x_position, t, dt);
        auto const mu =
            liquid_phase.property(MPL::PropertyType::viscosity)
                .template value<double>(variables, x_position, t, dt);

        // Set mechanical variables for the intrinsic permeability model
        // For stress dependent permeability.

        // For stress dependent permeability.
        variables[static_cast<int>(MPL::Variable::total_stress)]
            .emplace<SymmetricTensor>(
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                    sigma_total));

        variables[static_cast<int>(
            MaterialPropertyLib::Variable::equivalent_plastic_strain)] =
            ip_data.material_state_variables->getEquivalentPlasticStrain();

        auto const K_intrinsic = MPL::formEigenTensor<DisplacementDim>(
            medium->property(MPL::PropertyType::permeability)
                .value(variables, x_position, t, dt));

        GlobalDimMatrixType const Ki_over_mu = K_intrinsic / mu;
        GlobalDimMatrixType const rho_Ki_over_mu = rho_LR * Ki_over_mu;

        laplace_p.noalias() +=
            dNdx.transpose() * k_rel * rho_Ki_over_mu * dNdx * w;

        auto const beta_LR = 1 / rho_LR *
                             liquid_phase.property(MPL::PropertyType::density)
                                 .template dValue<double>(
                                     variables, MPL::Variable::phase_pressure,
                                     x_position, t, dt);

        const double alphaB_minus_phi = alpha - phi;
        double const a0 = alphaB_minus_phi * beta_SR;
        double const specific_storage_a_p = S_L * (phi * beta_LR + S_L * a0);
        double const specific_storage_a_S = phi - p_cap_ip * S_L * a0;

        // Note: d beta_LR/d p is omitted because it is a small value.
        double const dspecific_storage_a_p_dp_cap =
            dS_L_dp_cap * (phi * beta_LR + 2 * S_L * a0);
        double const dspecific_storage_a_S_dp_cap =
            -a0 * (S_L + p_cap_ip * dS_L_dp_cap);

        storage_p_a_p.noalias() +=
            N.transpose() * rho_LR * specific_storage_a_p * N * w;

        storage_p_a_S.noalias() -= N.transpose() * rho_LR *
                                   specific_storage_a_S * DeltaS_L_Deltap_cap *
                                   N * w;

        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() += N.transpose() * p_cap_dot_ip * rho_LR *
                          dspecific_storage_a_p_dp_cap * N * w;

        storage_p_a_S_Jpp.noalias() -=
            N.transpose() * rho_LR *
            ((S_L - S_L_prev) * dspecific_storage_a_S_dp_cap +
             specific_storage_a_S * dS_L_dp_cap) /
            dt * N * w;

        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() -= N.transpose() * rho_LR * dS_L_dp_cap * alpha *
                          identity2.transpose() * B * u_dot * N * w;

        double const dk_rel_dS_L =
            medium->property(MPL::PropertyType::relative_permeability)
                .template dValue<double>(variables,
                                         MPL::Variable::liquid_saturation,
                                         x_position, t, dt);
        GlobalDimVectorType const grad_p_cap = -dNdx * p_L;
        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() += dNdx.transpose() * rho_Ki_over_mu * grad_p_cap *
                          dk_rel_dS_L * dS_L_dp_cap * N * w;

        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() += dNdx.transpose() * rho_LR * rho_Ki_over_mu * b *
                          dk_rel_dS_L * dS_L_dp_cap * N * w;

        local_rhs.template segment<pressure_size>(pressure_index).noalias() +=
            dNdx.transpose() * rho_LR * k_rel * rho_Ki_over_mu * b * w;

        //
        // pressure equation, temperature part, thermal expansion.
        //
        {
            double const fluid_volumetric_thermal_expansion =
                phi * MPL::getLiquidThermalExpansivity(
                          liquid_phase, variables, rho_LR, x_position, t, dt);

            const double eff_thermal_expansion =
                alphaB_minus_phi *
                    Invariants::trace(solid_linear_thermal_expansivity_vector) +
                fluid_volumetric_thermal_expansion;

            M_pT.noalias() -=
                N.transpose() * S_L * rho_LR * eff_thermal_expansion * N * w;
        }

        //
        // temperature equation.
        //
        {
            auto const specific_heat_capacity_fluid =
                liquid_phase
                    .property(MaterialPropertyLib::specific_heat_capacity)
                    .template value<double>(variables, x_position, t, dt);

            auto const specific_heat_capacity_solid =
                solid_phase
                    .property(MaterialPropertyLib::PropertyType::
                                  specific_heat_capacity)
                    .template value<double>(variables, x_position, t, dt);

            M_TT.noalias() +=
                w *
                (rho_SR * specific_heat_capacity_solid * (1 - phi) +
                 (S_L * rho_LR * specific_heat_capacity_fluid) * phi) *
                N.transpose() * N;

            auto const thermal_conductivity =
                MaterialPropertyLib::formEigenTensor<DisplacementDim>(
                    medium
                        ->property(MaterialPropertyLib::PropertyType::
                                       thermal_conductivity)
                        .value(variables, x_position, t, dt));

            GlobalDimVectorType const velocity_L = GlobalDimVectorType(
                -Ki_over_mu * k_rel * (dNdx * p_L - rho_LR * b));

            K_TT.noalias() += (dNdx.transpose() * thermal_conductivity * dNdx +
                               N.transpose() * velocity_L.transpose() * dNdx *
                                   rho_LR * specific_heat_capacity_fluid) *
                              w;

            //
            // temperature equation, pressure part
            //
            K_Tp.noalias() -= rho_LR * specific_heat_capacity_fluid *
                              N.transpose() * (dNdx * T).transpose() * k_rel *
                              Ki_over_mu * dNdx * w;
            K_Tp.noalias() -= rho_LR * specific_heat_capacity_fluid *
                              N.transpose() * velocity_L.dot(dNdx * T) / k_rel *
                              dk_rel_dS_L * dS_L_dp_cap * N * w;
        }
    }

    if (process_data_.apply_mass_lumping)
    {
        storage_p_a_p = storage_p_a_p.colwise().sum().eval().asDiagonal();
        storage_p_a_S = storage_p_a_S.colwise().sum().eval().asDiagonal();
        storage_p_a_S_Jpp =
            storage_p_a_S_Jpp.colwise().sum().eval().asDiagonal();
    }

    //
    // -- Jacobian
    //
    // temperature equation.
    local_Jac
        .template block<temperature_size, temperature_size>(temperature_index,
                                                            temperature_index)
        .noalias() += M_TT / dt + K_TT;
    // temperature equation, pressure part
    local_Jac
        .template block<temperature_size, pressure_size>(temperature_index,
                                                         pressure_index)
        .noalias() += K_Tp;

    // pressure equation, pressure part.
    local_Jac
        .template block<pressure_size, pressure_size>(pressure_index,
                                                      pressure_index)
        .noalias() += laplace_p + storage_p_a_p / dt + storage_p_a_S_Jpp;

    // pressure equation, temperature part (contributed by thermal expansion).
    local_Jac
        .template block<pressure_size, temperature_size>(pressure_index,
                                                         temperature_index)
        .noalias() += M_pT / dt;

    // pressure equation, displacement part.
    local_Jac
        .template block<pressure_size, displacement_size>(pressure_index,
                                                          displacement_index)
        .noalias() = Kpu / dt;

    //
    // -- Residual
    //
    // temperature equation
    local_rhs.template segment<temperature_size>(temperature_index).noalias() -=
        M_TT * T_dot + K_TT * T;

    // pressure equation
    local_rhs.template segment<pressure_size>(pressure_index).noalias() -=
        laplace_p * p_L + (storage_p_a_p + storage_p_a_S) * p_L_dot +
        Kpu * u_dot + M_pT * T_dot;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunction, IntegrationMethod,
    DisplacementDim>::getSigma() const
{
    constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);

    return transposeInPlace<kelvin_vector_size>(
        [this](std::vector<double>& values)
        { return getIntPtSigma(0, {}, {}, values); });
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const&
ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement, ShapeFunction,
                                      IntegrationMethod, DisplacementDim>::
    getIntPtSigma(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
        ip_data_, &IpData::sigma_eff, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunction, IntegrationMethod,
    DisplacementDim>::getSwellingStress() const
{
    constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);

    return transposeInPlace<kelvin_vector_size>(
        [this](std::vector<double>& values)
        { return getIntPtSwellingStress(0, {}, {}, values); });
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const&
ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement, ShapeFunction,
                                      IntegrationMethod, DisplacementDim>::
    getIntPtSwellingStress(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    auto const n_integration_points = ip_data_.size();

    cache.clear();
    auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, kelvin_vector_size, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        auto const& sigma_sw = ip_data_[ip].sigma_sw;
        cache_mat.col(ip) =
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_sw);
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunction, IntegrationMethod,
    DisplacementDim>::getEpsilon() const
{
    constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);

    return transposeInPlace<kelvin_vector_size>(
        [this](std::vector<double>& values)
        { return getIntPtEpsilon(0, {}, {}, values); });
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const&
ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement, ShapeFunction,
                                      IntegrationMethod, DisplacementDim>::
    getIntPtEpsilon(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
        ip_data_, &IpData::eps, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const&
ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement, ShapeFunction,
                                      IntegrationMethod, DisplacementDim>::
    getIntPtDarcyVelocity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    unsigned const n_integration_points =
        integration_method_.getNumberOfPoints();

    cache.clear();
    auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, DisplacementDim, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        cache_matrix.col(ip).noalias() = ip_data_[ip].v_darcy;
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunction, IntegrationMethod,
    DisplacementDim>::getSaturation() const
{
    std::vector<double> result;
    getIntPtSaturation(0, {}, {}, result);
    return result;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const&
ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement, ShapeFunction,
                                      IntegrationMethod, DisplacementDim>::
    getIntPtSaturation(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(
        ip_data_, &IpData::saturation, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunction, IntegrationMethod,
    DisplacementDim>::getPorosity() const
{
    std::vector<double> result;
    getIntPtPorosity(0, {}, {}, result);
    return result;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const&
ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement, ShapeFunction,
                                      IntegrationMethod, DisplacementDim>::
    getIntPtPorosity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(ip_data_,
                                                     &IpData::porosity, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunction, IntegrationMethod,
    DisplacementDim>::getTransportPorosity() const
{
    std::vector<double> result;
    getIntPtTransportPorosity(0, {}, {}, result);
    return result;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const&
ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement, ShapeFunction,
                                      IntegrationMethod, DisplacementDim>::
    getIntPtTransportPorosity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(
        ip_data_, &IpData::transport_porosity, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const&
ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement, ShapeFunction,
                                      IntegrationMethod, DisplacementDim>::
    getIntPtDryDensitySolid(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(
        ip_data_, &IpData::dry_density_solid, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
void ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                           ShapeFunction, IntegrationMethod,
                                           DisplacementDim>::
    computeSecondaryVariableConcrete(double const t, double const dt,
                                     Eigen::VectorXd const& local_x,
                                     Eigen::VectorXd const& local_x_dot)
{
    auto const T =
        local_x.template segment<temperature_size>(temperature_index);
    auto const p_L = local_x.template segment<pressure_size>(pressure_index);
    auto const u =
        local_x.template segment<displacement_size>(displacement_index);

    auto const T_dot =
        local_x_dot.template segment<temperature_size>(temperature_index);
    auto const p_L_dot =
        local_x_dot.template segment<pressure_size>(pressure_index);
    auto const u_dot =
        local_x_dot.template segment<displacement_size>(displacement_index);

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    auto const& medium = process_data_.media_map->getMedium(element_.getID());
    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& solid_phase = medium->phase("Solid");
    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(element_.getID());

    unsigned const n_integration_points =
        integration_method_.getNumberOfPoints();

    double saturation_avg = 0;
    double porosity_avg = 0;

    using KV = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    KV sigma_avg = KV::Zero();

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto& ip_data = ip_data_[ip];
        // N is used for both p and T variables
        auto const& N = ip_data.N_p;
        auto const& N_u = ip_data.N_u;
        auto const& dNdx_u = ip_data.dNdx_u;

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                element_, N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, is_axially_symmetric_);

        double T_ip;
        NumLib::shapeFunctionInterpolate(T, N, T_ip);
        double T_dot_ip;
        NumLib::shapeFunctionInterpolate(T_dot, N, T_dot_ip);
        double const dT = T_dot_ip * dt;

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N, p_cap_ip);

        double p_cap_dot_ip;
        NumLib::shapeFunctionInterpolate(-p_L_dot, N, p_cap_dot_ip);

        variables[static_cast<int>(MPL::Variable::capillary_pressure)] =
            p_cap_ip;
        variables[static_cast<int>(MPL::Variable::phase_pressure)] = -p_cap_ip;

        variables[static_cast<int>(MPL::Variable::temperature)] = T_ip;

        auto& eps = ip_data.eps;
        eps.noalias() = B * u;
        auto& eps_m = ip_data.eps_m;
        auto& S_L = ip_data.saturation;
        auto const S_L_prev = ip_data.saturation_prev;
        S_L = medium->property(MPL::PropertyType::saturation)
                  .template value<double>(variables, x_position, t, dt);
        variables[static_cast<int>(MPL::Variable::liquid_saturation)] = S_L;
        variables_prev[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L_prev;

        auto const chi = [medium, x_position, t, dt](double const S_L)
        {
            MPL::VariableArray vs;
            vs.fill(std::numeric_limits<double>::quiet_NaN());
            vs[static_cast<int>(MPL::Variable::liquid_saturation)] = S_L;
            return medium->property(MPL::PropertyType::bishops_effective_stress)
                .template value<double>(vs, x_position, t, dt);
        };
        double const chi_S_L = chi(S_L);
        double const chi_S_L_prev = chi(S_L_prev);

        auto const alpha =
            medium->property(MPL::PropertyType::biot_coefficient)
                .template value<double>(variables, x_position, t, dt);

        double const T_ip_prev = T_ip - dT;
        auto const C_el = ip_data.computeElasticTangentStiffness(
            t, x_position, dt, T_ip_prev, T_ip);

        auto const beta_SR =
            (1 - alpha) /
            ip_data.solid_material.getBulkModulus(t, x_position, &C_el);
        variables[static_cast<int>(MPL::Variable::grain_compressibility)] =
            beta_SR;

        variables[static_cast<int>(MPL::Variable::effective_pore_pressure)] =
            -chi_S_L * p_cap_ip;
        variables_prev[static_cast<int>(
            MPL::Variable::effective_pore_pressure)] =
            -chi_S_L_prev * (p_cap_ip - p_cap_dot_ip * dt);

        // Set volumetric strain rate for the general case without swelling.
        variables[static_cast<int>(MPL::Variable::volumetric_strain)]
            .emplace<double>(Invariants::trace(eps));
        variables_prev[static_cast<int>(MPL::Variable::volumetric_strain)]
            .emplace<double>(Invariants::trace(B * (u - u_dot * dt)));

        auto& phi = ip_data.porosity;
        {  // Porosity update
            variables_prev[static_cast<int>(MPL::Variable::porosity)] =
                ip_data.porosity_prev;
            phi = medium->property(MPL::PropertyType::porosity)
                      .template value<double>(variables, variables_prev,
                                              x_position, t, dt);
            variables[static_cast<int>(MPL::Variable::porosity)] = phi;
        }

        // Swelling and possibly volumetric strain rate update.
        if (solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate))
        {
            auto& sigma_sw = ip_data.sigma_sw;
            auto const& sigma_sw_prev = ip_data.sigma_sw_prev;

            // If there is swelling, compute it. Update volumetric strain rate,
            // s.t. it corresponds to the mechanical part only.
            sigma_sw = sigma_sw_prev;

            using DimMatrix = Eigen::Matrix<double, 3, 3>;
            auto const sigma_sw_dot =
                MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
                    solid_phase
                        .property(MPL::PropertyType::swelling_stress_rate)
                        .template value<DimMatrix>(variables, variables_prev,
                                                   x_position, t, dt));
            sigma_sw += sigma_sw_dot * dt;

            // !!! Misusing volumetric strain for mechanical volumetric
            // strain just to update the transport porosity !!!
            std::get<double>(variables[static_cast<int>(
                MPL::Variable::volumetric_strain)]) +=
                identity2.transpose() * C_el.inverse() * sigma_sw;
            std::get<double>(variables_prev[static_cast<int>(
                MPL::Variable::volumetric_strain)]) +=
                identity2.transpose() * C_el.inverse() * sigma_sw_prev;
        }

        if (solid_phase.hasProperty(MPL::PropertyType::transport_porosity))
        {
            variables_prev[static_cast<int>(
                MPL::Variable::transport_porosity)] =
                ip_data.transport_porosity_prev;

            ip_data.transport_porosity =
                solid_phase.property(MPL::PropertyType::transport_porosity)
                    .template value<double>(variables, variables_prev,
                                            x_position, t, dt);
            variables[static_cast<int>(MPL::Variable::transport_porosity)] =
                ip_data.transport_porosity;
        }
        else
        {
            variables[static_cast<int>(MPL::Variable::transport_porosity)] =
                phi;
        }

        auto const mu =
            liquid_phase.property(MPL::PropertyType::viscosity)
                .template value<double>(variables, x_position, t, dt);
        auto const rho_LR =
            liquid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);

        // Set mechanical variables for the intrinsic permeability model
        // For stress dependent permeability.
        {
            auto const sigma_total =
                (ip_data.sigma_eff + alpha * chi_S_L * identity2 * p_cap_ip)
                    .eval();
            // For stress dependent permeability.
            variables[static_cast<int>(MPL::Variable::total_stress)]
                .emplace<SymmetricTensor>(
                    MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                        sigma_total));
        }

        variables[static_cast<int>(
            MaterialPropertyLib::Variable::equivalent_plastic_strain)] =
            ip_data.material_state_variables->getEquivalentPlasticStrain();
        auto const K_intrinsic = MPL::formEigenTensor<DisplacementDim>(
            medium->property(MPL::PropertyType::permeability)
                .value(variables, x_position, t, dt));

        double const k_rel =
            medium->property(MPL::PropertyType::relative_permeability)
                .template value<double>(variables, x_position, t, dt);

        GlobalDimMatrixType const K_over_mu = k_rel * K_intrinsic / mu;

        auto const& sigma_eff = ip_data.sigma_eff;
        double const p_FR = -chi_S_L * p_cap_ip;
        // p_SR
        variables[static_cast<int>(MPL::Variable::solid_grain_pressure)] =
            p_FR - Invariants::trace(sigma_eff) / (3 * (1 - phi));
        auto const rho_SR =
            solid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);
        ip_data.dry_density_solid = (1 - phi) * rho_SR;

        MathLib::KelvinVector::KelvinVectorType<
            DisplacementDim> const solid_linear_thermal_expansivity_vector =
            MPL::formKelvinVectorFromThermalExpansivity<DisplacementDim>(
                solid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_expansivity)
                    .value(variables, x_position, t, dt));

        MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const
            dthermal_strain = solid_linear_thermal_expansivity_vector * dT;

        auto& eps_prev = ip_data.eps_prev;
        auto& eps_m_prev = ip_data.eps_m_prev;
        eps_m.noalias() = eps_m_prev + eps - eps_prev - dthermal_strain;

        if (solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate))
        {
            eps_m.noalias() -=
                -C_el.inverse() * (ip_data.sigma_sw - ip_data.sigma_sw_prev);
        }

        variables[static_cast<int>(
                      MaterialPropertyLib::Variable::mechanical_strain)]
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m);

        ip_data.updateConstitutiveRelation(variables, t, x_position, dt,
                                           T_ip_prev);

        auto const& b = process_data_.specific_body_force;

        // Compute the velocity
        auto const& dNdx = ip_data.dNdx_p;
        ip_data.v_darcy.noalias() =
            -K_over_mu * dNdx * p_L + rho_LR * K_over_mu * b;

        saturation_avg += S_L;
        porosity_avg += phi;
        sigma_avg += sigma_eff;
    }
    saturation_avg /= n_integration_points;
    porosity_avg /= n_integration_points;
    sigma_avg /= n_integration_points;

    (*process_data_.element_saturation)[element_.getID()] = saturation_avg;
    (*process_data_.element_porosity)[element_.getID()] = porosity_avg;

    Eigen::Map<KV>(&(*process_data_.element_stresses)[element_.getID() *
                                                      KV::RowsAtCompileTime]) =
        MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_avg);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunction, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(element_, is_axially_symmetric_, p_L,
                         *process_data_.pressure_interpolated);
    NumLib::interpolateToHigherOrderNodes<
        ShapeFunction, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(element_, is_axially_symmetric_, T,
                         *process_data_.temperature_interpolated);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
unsigned ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunction, IntegrationMethod,
    DisplacementDim>::getNumberOfIntegrationPoints() const
{
    return integration_method_.getNumberOfPoints();
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
typename MaterialLib::Solids::MechanicsBase<
    DisplacementDim>::MaterialStateVariables const&
ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement, ShapeFunction,
                                      IntegrationMethod, DisplacementDim>::
    getMaterialStateVariablesAt(unsigned integration_point) const
{
    return *ip_data_[integration_point].material_state_variables;
}
}  // namespace ThermoRichardsMechanics
}  // namespace ProcessLib
