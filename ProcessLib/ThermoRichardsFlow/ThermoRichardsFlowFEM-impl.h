/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <cassert>

#include "HydrostaticElasticityModel.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEffectiveThermalConductivity.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/MPL/Utils/FormEigenVector.h"
#include "MaterialLib/MPL/Utils/GetLiquidThermalExpansivity.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "RigidElasticityModel.h"
#include "UniaxialElasticityModel.h"
#include "UserDefinedElasticityModel.h"

namespace ProcessLib
{
namespace ThermoRichardsFlow
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
ThermoRichardsFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    ThermoRichardsFlowLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        ThermoRichardsFlowProcessData& process_data)
    : _process_data(process_data),
      _integration_method(integration_order),
      _element(e),
      _is_axially_symmetric(is_axially_symmetric)
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);

    auto const shape_matrices =
        NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType, GlobalDim>(
            e, is_axially_symmetric, _integration_method);

    auto const& medium = *_process_data.media_map->getMedium(_element.getID());

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = shape_matrices[ip];
        x_position.setIntegrationPoint(ip);
        _ip_data.emplace_back();
        auto& ip_data = _ip_data[ip];
        _ip_data[ip].integration_weight =
            _integration_method.getWeightedPoint(ip).getWeight() *
            sm.integralMeasure * sm.detJ;

        ip_data.N = sm.N;
        ip_data.dNdx = sm.dNdx;

        // Initial porosity. Could be read from integration point data or mesh.
        ip_data.porosity =
            medium[MPL::porosity]
                .template initialValue<double>(
                    x_position,
                    std::numeric_limits<
                        double>::quiet_NaN() /* t independent */);
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
std::size_t ThermoRichardsFlowLocalAssembler<
    ShapeFunction, IntegrationMethod,
    GlobalDim>::setIPDataInitialConditions(std::string const& name,
                                           double const* values,
                                           int const integration_order)
{
    if (integration_order !=
        static_cast<int>(_integration_method.getIntegrationOrder()))
    {
        OGS_FATAL(
            "Setting integration point initial conditions; The integration "
            "order of the local assembler for element {:d} is different "
            "from the integration order in the initial condition.",
            _element.getID());
    }

    if (name == "saturation_ip")
    {
        return ProcessLib::setIntegrationPointScalarData(values, _ip_data,
                                                         &IpData::saturation);
    }
    if (name == "porosity_ip")
    {
        return ProcessLib::setIntegrationPointScalarData(values, _ip_data,
                                                         &IpData::porosity);
    }
    return 0;
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void ThermoRichardsFlowLocalAssembler<ShapeFunction, IntegrationMethod,
                                      GlobalDim>::
    setInitialConditionsConcrete(std::vector<double> const& local_x,
                                 double const t,
                                 bool const /*use_monolithic_scheme*/,
                                 int const /*process_id*/)
{
    assert(local_x.size() == temperature_size + pressure_size);

    auto p_L = Eigen::Map<
        typename ShapeMatricesType::template VectorType<pressure_size> const>(
        local_x.data() + pressure_index, pressure_size);

    auto const& medium = *_process_data.media_map->getMedium(_element.getID());
    MPL::VariableArray variables;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto const& N = _ip_data[ip].N;

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N, p_cap_ip);

        variables[static_cast<int>(MPL::Variable::capillary_pressure)] =
            p_cap_ip;
        variables[static_cast<int>(MPL::Variable::phase_pressure)] = -p_cap_ip;

        // Note: temperature dependent saturation model is not considered so
        // far.
        _ip_data[ip].saturation_prev =
            medium[MPL::PropertyType::saturation]
                .template value<double>(
                    variables, x_position, t,
                    std::numeric_limits<double>::quiet_NaN());
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void ThermoRichardsFlowLocalAssembler<
    ShapeFunction, IntegrationMethod,
    GlobalDim>::assembleWithJacobian(double const t, double const dt,
                                     std::vector<double> const& local_x,
                                     std::vector<double> const& local_xdot,
                                     const double /*dxdot_dx*/,
                                     const double /*dx_dx*/,
                                     std::vector<double>& /*local_M_data*/,
                                     std::vector<double>& /*local_K_data*/,
                                     std::vector<double>& local_rhs_data,
                                     std::vector<double>& local_Jac_data)
{
    auto const local_matrix_dim = pressure_size + temperature_size;
    assert(local_x.size() == local_matrix_dim);

    auto const T = Eigen::Map<typename ShapeMatricesType::template VectorType<
        temperature_size> const>(local_x.data() + temperature_index,
                                 temperature_size);
    auto const p_L = Eigen::Map<
        typename ShapeMatricesType::template VectorType<pressure_size> const>(
        local_x.data() + pressure_index, pressure_size);

    auto const T_dot =
        Eigen::Map<typename ShapeMatricesType::template VectorType<
            temperature_size> const>(local_xdot.data() + temperature_index,
                                     temperature_size);
    auto const p_L_dot = Eigen::Map<
        typename ShapeMatricesType::template VectorType<pressure_size> const>(
        local_xdot.data() + pressure_index, pressure_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesType::template MatrixType<local_matrix_dim,
                                                        local_matrix_dim>>(
        local_Jac_data, local_matrix_dim, local_matrix_dim);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesType::template VectorType<local_matrix_dim>>(
        local_rhs_data, local_matrix_dim);

    typename ShapeMatricesType::NodalMatrixType M_TT =
        ShapeMatricesType::NodalMatrixType::Zero(temperature_size,
                                                 temperature_size);
    typename ShapeMatricesType::NodalMatrixType M_Tp =
        ShapeMatricesType::NodalMatrixType::Zero(temperature_size,
                                                 pressure_size);
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

    auto const& medium = *_process_data.media_map->getMedium(_element.getID());
    auto const& liquid_phase = medium.phase("AqueousLiquid");
    auto const& solid_phase = medium.phase("Solid");
    MPL::Phase const* gas_phase = medium.hasPhase("Gas") ?
        &medium.phase("Gas") : nullptr;
    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;

        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        double T_ip;
        NumLib::shapeFunctionInterpolate(T, N, T_ip);
        double T_dot_ip;
        NumLib::shapeFunctionInterpolate(T_dot, N, T_dot_ip);

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N, p_cap_ip);

        double p_cap_dot_ip;
        NumLib::shapeFunctionInterpolate(-p_L_dot, N, p_cap_dot_ip);

        variables[static_cast<int>(MPL::Variable::capillary_pressure)] =
            p_cap_ip;
        variables[static_cast<int>(MPL::Variable::phase_pressure)] = -p_cap_ip;
        variables[static_cast<int>(MPL::Variable::temperature)] = T_ip;

        auto& S_L = _ip_data[ip].saturation;
        auto const S_L_prev = _ip_data[ip].saturation_prev;
        auto const alpha =
            medium[MPL::PropertyType::biot_coefficient]
                .template value<double>(variables, x_position, t, dt);

        auto& solid_elasticity = *_process_data.simplified_elasticity;
        // TODO (buchwaldj)
        // is bulk_modulus good name for bulk modulus of solid skeleton?
        auto const beta_S =
            solid_elasticity.bulkCompressibilityFromYoungsModulus(
                solid_phase, variables, x_position, t, dt);
        auto const beta_SR = (1 - alpha) * beta_S;
        variables[static_cast<int>(MPL::Variable::grain_compressibility)] =
            beta_SR;

        auto const rho_LR =
            liquid_phase[MPL::PropertyType::density]
                .template value<double>(variables, x_position, t, dt);
        auto const& b = _process_data.specific_body_force;

        double const drho_LR_dp =
            liquid_phase[MPL::PropertyType::density]
                .template dValue<double>(variables,
                                         MPL::Variable::phase_pressure,
                                         x_position, t, dt);
        auto const beta_LR = drho_LR_dp / rho_LR;

        S_L = medium[MPL::PropertyType::saturation]
                  .template value<double>(variables, x_position, t, dt);
        variables[static_cast<int>(MPL::Variable::liquid_saturation)] = S_L;
        variables_prev[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L_prev;

        // tangent derivative for Jacobian
        double const dS_L_dp_cap = medium[MPL::PropertyType::saturation]
                .template dValue<double>(variables,
                                         MPL::Variable::capillary_pressure,
                                         x_position, t, dt);
        // secant derivative from time discretization for storage
        // use tangent, if secant is not available
        double const DeltaS_L_Deltap_cap =
            (p_cap_dot_ip == 0) ? dS_L_dp_cap
                                : (S_L - S_L_prev) / (dt * p_cap_dot_ip);

        auto chi_S_L = S_L;
        auto chi_S_L_prev = S_L_prev;
        auto dchi_dS_L = 1.0;
        if (medium.hasProperty(MPL::PropertyType::bishops_effective_stress))
        {
            auto const chi = [&medium, x_position, t, dt](double const S_L) {
                MPL::VariableArray variables;
                variables[static_cast<int>(MPL::Variable::liquid_saturation)] =
                    S_L;
                return medium[MPL::PropertyType::bishops_effective_stress]
                    .template value<double>(variables, x_position, t, dt);
            };
            chi_S_L = chi(S_L);
            chi_S_L_prev = chi(S_L_prev);

            dchi_dS_L = medium[MPL::PropertyType::bishops_effective_stress]
                .template dValue<double>(variables,
                                         MPL::Variable::liquid_saturation,
                                         x_position, t, dt);
        }
        // TODO (buchwaldj)
        // should solid_grain_pressure or effective_pore_pressure remain?
        // double const p_FR = -chi_S_L * p_cap_ip;
        // variables[static_cast<int>(MPL::Variable::solid_grain_pressure)] =
        // p_FR;

        variables[static_cast<int>(MPL::Variable::effective_pore_pressure)] =
            -chi_S_L * p_cap_ip;
        variables_prev[static_cast<int>(
            MPL::Variable::effective_pore_pressure)] =
            -chi_S_L_prev * (p_cap_ip - p_cap_dot_ip * dt);

        auto& phi = _ip_data[ip].porosity;
        {  // Porosity update

            variables_prev[static_cast<int>(MPL::Variable::porosity)] =
                _ip_data[ip].porosity_prev;
            phi = medium[MPL::PropertyType::porosity]
                      .template value<double>(variables, variables_prev,
                                              x_position, t, dt);
            variables[static_cast<int>(MPL::Variable::porosity)] = phi;
        }

        if (alpha < phi)
        {
            OGS_FATAL(
                "ThermoRichardsFlow: Biot-coefficient {} is smaller than "
                "porosity {} in element/integration point {}/{}.",
                alpha, phi, _element.getID(), ip);
        }

        double const k_rel =
            medium[MPL::PropertyType::relative_permeability]
                .template value<double>(variables, x_position, t, dt);
        auto const mu =
            liquid_phase[MPL::PropertyType::viscosity]
                .template value<double>(variables, x_position, t, dt);

        auto const K_intrinsic = MPL::formEigenTensor<GlobalDim>(
            medium[MPL::PropertyType::permeability]
                .value(variables, x_position, t, dt));

        GlobalDimMatrixType const Ki_over_mu = K_intrinsic / mu;
        GlobalDimMatrixType const rho_Ki_over_mu = rho_LR * Ki_over_mu;

        // Consider anisotropic thermal expansion.
        // Read in 3x3 tensor. 2D case also requires expansion coeff. for z-
        // component.
        Eigen::Matrix<double, 3,
                      3> const solid_linear_thermal_expansion_coefficient =
            MaterialPropertyLib::formEigenTensor<3>(
                solid_phase[
                        MaterialPropertyLib::PropertyType::thermal_expansivity]
                    .value(variables, x_position, t, dt));

        auto const rho_SR =
            solid_phase[MPL::PropertyType::density]
                .template value<double>(variables, x_position, t, dt);

        //
        // pressure equation, pressure part.
        //
        laplace_p.noalias() +=
            dNdx.transpose() * k_rel * rho_Ki_over_mu * dNdx * w;

        const double alphaB_minus_phi = alpha - phi;
        double const a0 = alphaB_minus_phi * beta_SR;
        double const specific_storage_a_p =
            S_L * (phi * beta_LR + S_L * a0 +
                   chi_S_L * alpha * alpha *
                       solid_elasticity.storageContribution(
                           solid_phase, variables, x_position, t, dt));
        double const specific_storage_a_S = phi - p_cap_ip * S_L * a0;

        double const dspecific_storage_a_p_dp_cap =
            dS_L_dp_cap * (phi * beta_LR + 2 * S_L * a0 +
                           alpha * alpha *
                               solid_elasticity.storageContribution(
                                   solid_phase, variables, x_position, t, dt) *
                               (chi_S_L + dchi_dS_L * S_L));
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

        double const dk_rel_dS_L =
            medium[MPL::PropertyType::relative_permeability]
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
        // pressure equation, temperature part.
        //
        double const fluid_volumetric_thermal_expansion_coefficient =
            MPL::getLiquidThermalExpansivity(liquid_phase, variables, rho_LR,
                                             x_position, t, dt);
        const double eff_thermal_expansion =
            S_L * (alphaB_minus_phi *
                       solid_linear_thermal_expansion_coefficient.trace() +
                   phi * fluid_volumetric_thermal_expansion_coefficient +
                   alpha * solid_elasticity.thermalExpansivityContribution(
                               solid_linear_thermal_expansion_coefficient,
                               solid_phase, variables, x_position, t, dt));
        M_pT.noalias() -=
            N.transpose() * rho_LR * eff_thermal_expansion * N * w;

        //
        // temperature equation.
        //
        {
            auto const specific_heat_capacity_fluid =
                liquid_phase[MaterialPropertyLib::specific_heat_capacity]
                    .template value<double>(variables, x_position, t, dt);

            auto const specific_heat_capacity_solid =
                solid_phase
                    [MaterialPropertyLib::PropertyType::
                                  specific_heat_capacity]
                    .template value<double>(variables, x_position, t, dt);

            M_TT.noalias() +=
                w *
                (rho_SR * specific_heat_capacity_solid * (1 - phi) +
                 (S_L * rho_LR * specific_heat_capacity_fluid) * phi) *
                N.transpose() * N;

            auto const thermal_conductivity =
                MaterialPropertyLib::formEigenTensor<GlobalDim>(
                            medium[MaterialPropertyLib::PropertyType::
                                               thermal_conductivity]
                                .value(variables, x_position, t, dt));

            GlobalDimVectorType const velocity_L =
                GlobalDimVectorType(-Ki_over_mu * (dNdx * p_L - rho_LR * b));

            K_TT.noalias() +=
                (dNdx.transpose() * thermal_conductivity * dNdx +
                 N.transpose() * velocity_L.transpose() * dNdx * rho_LR *
                     specific_heat_capacity_fluid) *
                w;

            //
            // temperature equation, pressure part
            //
            K_Tp.noalias() -= rho_LR * specific_heat_capacity_fluid *
                              N.transpose() * (dNdx * T).transpose() *
                              Ki_over_mu * dNdx * w;
        }
        if (liquid_phase.hasProperty(MPL::PropertyType::vapour_diffusion) &&
            S_L < 1.0)
        {
            variables[static_cast<int>(MPL::Variable::density)] = rho_LR;

            double const rho_wv =
                liquid_phase[MaterialPropertyLib::vapour_density]
                    .template value<double>(variables, x_position, t, dt);

            double const drho_wv_dT =
                liquid_phase[MaterialPropertyLib::vapour_density]
                    .template dValue<double>(variables,
                                             MPL::Variable::temperature,
                                             x_position, t, dt);
            double const drho_wv_dp =
                liquid_phase[MaterialPropertyLib::vapour_density]
                    .template dValue<double>(variables,
                                             MPL::Variable::phase_pressure,
                                             x_position, t, dt);
            auto const f_Tv =
                liquid_phase[
                    MPL::PropertyType::thermal_diffusion_enhancement_factor]
                    .template value<double>(variables, x_position, t, dt);

            variables[static_cast<int>(MPL::Variable::porosity)] = phi;
            double const D_v =
                liquid_phase[MPL::PropertyType::vapour_diffusion]
                    .template value<double>(variables, x_position, t, dt);

            double const f_Tv_D_Tv = f_Tv * D_v * drho_wv_dT;
            double const D_pv = D_v * drho_wv_dp;

            if (gas_phase &&
                gas_phase->hasProperty(MPL::PropertyType::heat_capacity))
            {
                GlobalDimVectorType const grad_T = dNdx * T;
                // Vapour velocity
                GlobalDimVectorType const q_v =
                    -(f_Tv_D_Tv * grad_T - D_pv * grad_p_cap) / rho_LR;
                double const specific_heat_capacity_vapour =
                    gas_phase->property(MaterialPropertyLib::PropertyType::
                                       specific_heat_capacity)
                        .template value<double>(variables, x_position, t, dt);

                M_TT.noalias() +=
                    w *
                    (rho_wv * specific_heat_capacity_vapour * (1 - S_L) * phi) *
                    N.transpose() * N;

                K_TT.noalias() += N.transpose() * q_v.transpose() * dNdx *
                                  rho_wv * specific_heat_capacity_vapour * w;
            }

            double const storage_coefficient_by_water_vapor =
                phi * (rho_wv * dS_L_dp_cap + (1 - S_L) * drho_wv_dp);

            storage_p_a_p.noalias() +=
                N.transpose() * storage_coefficient_by_water_vapor * N * w;

            double const vapor_expansion_factor = phi * (1 - S_L) * drho_wv_dT;
            M_pT.noalias() += N.transpose() * vapor_expansion_factor * N * w;

            local_Jac
                .template block<pressure_size, temperature_size>(
                    pressure_index, temperature_index)
                .noalias() += dNdx.transpose() * f_Tv_D_Tv * dNdx * w;

            local_rhs.template segment<pressure_size>(pressure_index)
                .noalias() -= f_Tv_D_Tv * dNdx.transpose() * (dNdx * T) * w;

            laplace_p.noalias() += dNdx.transpose() * D_pv * dNdx * w;

            //
            // Latent heat term
            //
            if (liquid_phase.hasProperty(MPL::PropertyType::latent_heat))
            {
                double const factor = phi * (1 - S_L) / rho_LR;
                // The volumetric latent heat of vaporization of liquid water
                double const L0 =
                    liquid_phase[MPL::PropertyType::latent_heat]
                        .template value<double>(variables, x_position, t, dt) *
                    rho_LR;

                double const drho_LR_dT =
                    liquid_phase[MPL::PropertyType::density]
                        .template dValue<double>(variables,
                                                 MPL::Variable::temperature,
                                                 x_position, t, dt);

                double const rho_wv_over_rho_L = rho_wv / rho_LR;
                M_TT.noalias() +=
                    factor * L0 *
                    (drho_wv_dT - rho_wv_over_rho_L * drho_LR_dT) *
                    N.transpose() * N * w;

                M_Tp.noalias() +=
                    (factor * L0 *
                         (drho_wv_dp - rho_wv_over_rho_L * drho_LR_dp) +
                     L0 * phi * rho_wv_over_rho_L * dS_L_dp_cap) *
                    N.transpose() * N * w;

                // temperature equation, temperature part
                K_TT.noalias() +=
                    L0 * f_Tv_D_Tv * dNdx.transpose() * dNdx * w / rho_LR;
                // temperature equation, pressure part
                K_Tp.noalias() +=
                    L0 * D_pv * dNdx.transpose() * dNdx * w / rho_LR;
            }
        }
    }

    if (_process_data.apply_mass_lumping)
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

    //
    // -- Residual
    //
    // temperature equation
    local_rhs.template segment<temperature_size>(temperature_index).noalias() -=
        M_TT * T_dot + K_TT * T;

    // pressure equation
    local_rhs.template segment<pressure_size>(pressure_index).noalias() -=
        laplace_p * p_L + (storage_p_a_p + storage_p_a_S) * p_L_dot +
        M_pT * T_dot;
    if (liquid_phase.hasProperty(MPL::PropertyType::vapour_diffusion) &&
        liquid_phase.hasProperty(MPL::PropertyType::latent_heat))
    {
        // Jacobian: temperature equation, pressure part
        local_Jac
            .template block<temperature_size, pressure_size>(temperature_index,
                                                             pressure_index)
            .noalias() += M_Tp / dt;
        // RHS: temperature part
        local_rhs.template segment<temperature_size>(temperature_index)
            .noalias() -= M_Tp * p_L_dot;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void ThermoRichardsFlowLocalAssembler<
    ShapeFunction, IntegrationMethod,
    GlobalDim>::assemble(double const t, double const dt,
                                     std::vector<double> const& local_x,
                                     std::vector<double> const& local_xdot,
                                     std::vector<double>& local_M_data,
                                     std::vector<double>& local_K_data,
                                     std::vector<double>& local_rhs_data)
{
    auto const local_matrix_dim = pressure_size + temperature_size;
    assert(local_x.size() == local_matrix_dim);

    auto const T = Eigen::Map<typename ShapeMatricesType::template VectorType<
        temperature_size> const>(local_x.data() + temperature_index,
                                 temperature_size);
    auto const p_L = Eigen::Map<
        typename ShapeMatricesType::template VectorType<pressure_size> const>(
        local_x.data() + pressure_index, pressure_size);

    auto const T_dot =
        Eigen::Map<typename ShapeMatricesType::template VectorType<
            temperature_size> const>(local_xdot.data() + temperature_index,
                                     temperature_size);
    auto const p_L_dot = Eigen::Map<
        typename ShapeMatricesType::template VectorType<pressure_size> const>(
        local_xdot.data() + pressure_index, pressure_size);

    auto local_K = MathLib::createZeroedMatrix<
        typename ShapeMatricesType::template MatrixType<
            local_matrix_dim, local_matrix_dim>>(
        local_K_data, local_matrix_dim, local_matrix_dim);

    auto local_M = MathLib::createZeroedMatrix<
        typename ShapeMatricesType::template MatrixType<
            local_matrix_dim, local_matrix_dim>>(
        local_M_data, local_matrix_dim, local_matrix_dim);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesType::template VectorType<local_matrix_dim>>(
        local_rhs_data, local_matrix_dim);

    auto const& medium = *_process_data.media_map->getMedium(_element.getID());
    auto const& liquid_phase = medium.phase("AqueousLiquid");
    auto const& solid_phase = medium.phase("Solid");
    MPL::Phase const* gas_phase = medium.hasPhase("Gas") ?
        &medium.phase("Gas") : nullptr;
    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;

        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        double T_ip;
        NumLib::shapeFunctionInterpolate(T, N, T_ip);
        double T_dot_ip;
        NumLib::shapeFunctionInterpolate(T_dot, N, T_dot_ip);

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N, p_cap_ip);

        double p_cap_dot_ip;
        NumLib::shapeFunctionInterpolate(-p_L_dot, N, p_cap_dot_ip);

        variables[static_cast<int>(MPL::Variable::capillary_pressure)] =
            p_cap_ip;
        variables[static_cast<int>(MPL::Variable::phase_pressure)] = -p_cap_ip;
        variables[static_cast<int>(MPL::Variable::temperature)] = T_ip;

        auto& S_L = _ip_data[ip].saturation;
        auto const S_L_prev = _ip_data[ip].saturation_prev;
        auto const alpha =
            medium[MPL::PropertyType::biot_coefficient]
                .template value<double>(variables, x_position, t, dt);

        auto& solid_elasticity = *_process_data.simplified_elasticity;
        // TODO (buchwaldj)
        // is bulk_modulus good name for bulk modulus of solid skeleton?
        auto const beta_S =
            solid_elasticity.bulkCompressibilityFromYoungsModulus(
                solid_phase, variables, x_position, t, dt);
        auto const beta_SR = (1 - alpha) * beta_S;
        variables[static_cast<int>(MPL::Variable::grain_compressibility)] =
            beta_SR;

        auto const rho_LR =
            liquid_phase[MPL::PropertyType::density]
                .template value<double>(variables, x_position, t, dt);
        auto const& b = _process_data.specific_body_force;

        double const drho_LR_dp =
            liquid_phase[MPL::PropertyType::density]
                .template dValue<double>(variables,
                                         MPL::Variable::phase_pressure,
                                         x_position, t, dt);
        auto const beta_LR = drho_LR_dp / rho_LR;

        S_L = medium[MPL::PropertyType::saturation]
                  .template value<double>(variables, x_position, t, dt);
        variables[static_cast<int>(MPL::Variable::liquid_saturation)] = S_L;
        variables_prev[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L_prev;

        // tangent derivative for Jacobian
        double const dS_L_dp_cap = medium[MPL::PropertyType::saturation]
                .template dValue<double>(variables,
                                         MPL::Variable::capillary_pressure,
                                         x_position, t, dt);
        // secant derivative from time discretization for storage
        // use tangent, if secant is not available
        double const DeltaS_L_Deltap_cap =
            (p_cap_dot_ip == 0) ? dS_L_dp_cap
                                : (S_L - S_L_prev) / (dt * p_cap_dot_ip);

        auto chi_S_L = S_L;
        auto chi_S_L_prev = S_L_prev;
        if (medium.hasProperty(MPL::PropertyType::bishops_effective_stress))
        {
            auto const chi = [&medium, x_position, t, dt](double const S_L) {
                MPL::VariableArray variables;
                variables[static_cast<int>(MPL::Variable::liquid_saturation)] =
                    S_L;
                return medium[MPL::PropertyType::bishops_effective_stress]
                    .template value<double>(variables, x_position, t, dt);
            };
            chi_S_L = chi(S_L);
            chi_S_L_prev = chi(S_L_prev);

        }
        // TODO (buchwaldj)
        // should solid_grain_pressure or effective_pore_pressure remain?
        // double const p_FR = -chi_S_L * p_cap_ip;
        // variables[static_cast<int>(MPL::Variable::solid_grain_pressure)] =
        // p_FR;

        variables[static_cast<int>(MPL::Variable::effective_pore_pressure)] =
            -chi_S_L * p_cap_ip;
        variables_prev[static_cast<int>(
            MPL::Variable::effective_pore_pressure)] =
            -chi_S_L_prev * (p_cap_ip - p_cap_dot_ip * dt);

        auto& phi = _ip_data[ip].porosity;
        {  // Porosity update

            variables_prev[static_cast<int>(MPL::Variable::porosity)] =
                _ip_data[ip].porosity_prev;
            phi = medium[MPL::PropertyType::porosity]
                      .template value<double>(variables, variables_prev,
                                              x_position, t, dt);
            variables[static_cast<int>(MPL::Variable::porosity)] = phi;
        }

        if (alpha < phi)
        {
            OGS_FATAL(
                "ThermoRichardsFlow: Biot-coefficient {} is smaller than "
                "porosity {} in element/integration point {}/{}.",
                alpha, phi, _element.getID(), ip);
        }

        double const k_rel =
            medium[MPL::PropertyType::relative_permeability]
                .template value<double>(variables, x_position, t, dt);
        auto const mu =
            liquid_phase[MPL::PropertyType::viscosity]
                .template value<double>(variables, x_position, t, dt);

        auto const K_intrinsic = MPL::formEigenTensor<GlobalDim>(
            medium[MPL::PropertyType::permeability]
                .value(variables, x_position, t, dt));

        GlobalDimMatrixType const Ki_over_mu = K_intrinsic / mu;
        GlobalDimMatrixType const rho_Ki_over_mu = rho_LR * Ki_over_mu;

        // Consider anisotropic thermal expansion.
        // Read in 3x3 tensor. 2D case also requires expansion coeff. for z-
        // component.
        Eigen::Matrix<double, 3,
                      3> const solid_linear_thermal_expansion_coefficient =
            MaterialPropertyLib::formEigenTensor<3>(
                solid_phase[
                        MaterialPropertyLib::PropertyType::thermal_expansivity]
                    .value(variables, x_position, t, dt));

        auto const rho_SR =
            solid_phase[MPL::PropertyType::density]
                .template value<double>(variables, x_position, t, dt);

        //
        // pressure equation, pressure part.
        //
        local_K.template block<pressure_size, pressure_size>(pressure_index,
                                                       pressure_index)
            .noalias() += dNdx.transpose() * k_rel * rho_Ki_over_mu * dNdx * w;

        const double alphaB_minus_phi = alpha - phi;
        double const a0 = alphaB_minus_phi * beta_SR;
        double const specific_storage_a_p =
            S_L * (phi * beta_LR + S_L * a0 +
                   chi_S_L * alpha * alpha *
                       solid_elasticity.storageContribution(
                           solid_phase, variables, x_position, t, dt));
        double const specific_storage_a_S = phi - p_cap_ip * S_L * a0;


        local_M.template block<pressure_size, pressure_size>(pressure_index,
                                                       pressure_index)
            .noalias() += N.transpose() * rho_LR * (specific_storage_a_p -
                    specific_storage_a_S * DeltaS_L_Deltap_cap) * N * w;


        local_rhs.template segment<pressure_size>(pressure_index).noalias() +=
            dNdx.transpose() * rho_LR * k_rel * rho_Ki_over_mu * b * w;

        //
        // pressure equation, temperature part.
        //
        double const fluid_volumetric_thermal_expansion_coefficient =
            MPL::getLiquidThermalExpansivity(liquid_phase, variables, rho_LR,
                                             x_position, t, dt);
        const double eff_thermal_expansion =
            S_L * (alphaB_minus_phi *
                       solid_linear_thermal_expansion_coefficient.trace() +
                   phi * fluid_volumetric_thermal_expansion_coefficient +
                   alpha * solid_elasticity.thermalExpansivityContribution(
                               solid_linear_thermal_expansion_coefficient,
                               solid_phase, variables, x_position, t, dt));

        local_M.template block<pressure_size, temperature_size>(pressure_index,
                                                       temperature_index)
            .noalias() -=
            N.transpose() * rho_LR * eff_thermal_expansion * N * w;

        //
        // temperature equation.
        //
        {
            auto const specific_heat_capacity_fluid =
                liquid_phase[MaterialPropertyLib::specific_heat_capacity]
                    .template value<double>(variables, x_position, t, dt);

            auto const specific_heat_capacity_solid =
                solid_phase[MaterialPropertyLib::PropertyType::
                                  specific_heat_capacity]
                    .template value<double>(variables, x_position, t, dt);

            local_M.template block<temperature_size, temperature_size>(temperature_index,
                                                       temperature_index)
            .noalias() +=
                w *
                (rho_SR * specific_heat_capacity_solid * (1 - phi) +
                 (S_L * rho_LR * specific_heat_capacity_fluid) * phi) *
                N.transpose() * N;

            auto const thermal_conductivity =
                        MaterialPropertyLib::formEigenTensor<GlobalDim>(
                            medium[MaterialPropertyLib::PropertyType::
                                               thermal_conductivity]
                                .value(variables, x_position, t, dt));

            GlobalDimVectorType const velocity_L =
                GlobalDimVectorType(-Ki_over_mu * (dNdx * p_L - rho_LR * b));

            local_K.template block<temperature_size, temperature_size>(temperature_index,
                                                       temperature_index)
            .noalias() +=
                (dNdx.transpose() * thermal_conductivity * dNdx +
                 N.transpose() * velocity_L.transpose() * dNdx * rho_LR *
                     specific_heat_capacity_fluid) *
                w;

            //
            // temperature equation, pressure part
            //
            local_K.template block<temperature_size, pressure_size>(temperature_index,
                                                       pressure_index)
                .noalias() -=
                rho_LR * specific_heat_capacity_fluid *
                              N.transpose() * (dNdx * T).transpose() *
                              Ki_over_mu * dNdx * w;
        }
        if (liquid_phase.hasProperty(MPL::PropertyType::vapour_diffusion) &&
            S_L < 1.0)
        {
            variables[static_cast<int>(MPL::Variable::density)] = rho_LR;

            double const rho_wv =
                liquid_phase[MaterialPropertyLib::vapour_density]
                    .template value<double>(variables, x_position, t, dt);

            double const drho_wv_dT =
                liquid_phase[MaterialPropertyLib::vapour_density]
                    .template dValue<double>(variables,
                                             MPL::Variable::temperature,
                                             x_position, t, dt);
            double const drho_wv_dp =
                liquid_phase[MaterialPropertyLib::vapour_density]
                    .template dValue<double>(variables,
                                             MPL::Variable::phase_pressure,
                                             x_position, t, dt);
            auto const f_Tv =
                liquid_phase[
                        MPL::PropertyType::thermal_diffusion_enhancement_factor]
                    .template value<double>(variables, x_position, t, dt);

            variables[static_cast<int>(MPL::Variable::porosity)] = phi;
            double const D_v =
                liquid_phase[MPL::PropertyType::vapour_diffusion]
                    .template value<double>(variables, x_position, t, dt);

            double const f_Tv_D_Tv = f_Tv * D_v * drho_wv_dT;
            double const D_pv = D_v * drho_wv_dp;

            if (gas_phase &&
                gas_phase->hasProperty(MPL::PropertyType::heat_capacity))
            {
                GlobalDimVectorType const grad_T = dNdx * T;
                GlobalDimVectorType const grad_p_cap = -dNdx * p_L;
                // Vapour velocity
                GlobalDimVectorType const q_v =
                    -(f_Tv_D_Tv * grad_T - D_pv * grad_p_cap) / rho_LR;
                double const specific_heat_capacity_vapour =
                    gas_phase->property(
                        MaterialPropertyLib::PropertyType::
                                       specific_heat_capacity)
                        .template value<double>(variables, x_position, t, dt);

                local_M.template block<temperature_size, temperature_size>(temperature_index,
                                                       temperature_index)
                .noalias() +=
                    w *
                    (rho_wv * specific_heat_capacity_vapour * (1 - S_L) * phi) *
                    N.transpose() * N;

                local_K.template block<temperature_size, temperature_size>(temperature_index,
                                                       temperature_index)
                .noalias() += N.transpose() * q_v.transpose() * dNdx *
                                  rho_wv * specific_heat_capacity_vapour * w;
            }

            double const storage_coefficient_by_water_vapor =
                phi * (rho_wv * dS_L_dp_cap + (1 - S_L) * drho_wv_dp);
            local_M.template block<pressure_size, pressure_size>(pressure_index,
                                                       pressure_index)
            .noalias() +=
                N.transpose() * storage_coefficient_by_water_vapor * N * w;

            double const vapor_expansion_factor = phi * (1 - S_L) * drho_wv_dT;
            local_M.template block<pressure_size, temperature_size>(pressure_index,
                                                       temperature_index)
            .noalias() += N.transpose() * vapor_expansion_factor * N * w;

            local_rhs.template segment<pressure_size>(pressure_index)
                .noalias() -= f_Tv_D_Tv * dNdx.transpose() * (dNdx * T) * w;

            local_K.template block<pressure_size, pressure_size>(pressure_index,
                                                       pressure_index)
            .noalias() += dNdx.transpose() * D_pv * dNdx * w;

            //
            // Latent heat term
            //
            if (liquid_phase.hasProperty(MPL::PropertyType::latent_heat))
            {
                double const factor = phi * (1 - S_L) / rho_LR;
                // The volumetric latent heat of vaporization of liquid water
                double const L0 =
                    liquid_phase[MPL::PropertyType::latent_heat]
                        .template value<double>(variables, x_position, t, dt) *
                    rho_LR;

                double const drho_LR_dT =
                    liquid_phase[MPL::PropertyType::density]
                        .template dValue<double>(variables,
                                                 MPL::Variable::temperature,
                                                 x_position, t, dt);

                double const rho_wv_over_rho_L = rho_wv / rho_LR;
                local_M.template block<temperature_size, temperature_size>(temperature_index,
                                                       temperature_index)
                .noalias() +=
                    factor * L0 *
                    (drho_wv_dT - rho_wv_over_rho_L * drho_LR_dT) *
                    N.transpose() * N * w;

                local_M.template block<temperature_size, pressure_size>(temperature_index,
                                                       pressure_index)
                .noalias() +=
                    (factor * L0 *
                         (drho_wv_dp - rho_wv_over_rho_L * drho_LR_dp) +
                     L0 * phi * rho_wv_over_rho_L * dS_L_dp_cap) *
                    N.transpose() * N * w;

                // temperature equation, temperature part
                local_K.template block<temperature_size, temperature_size>(temperature_index,
                                                       temperature_index)
                .noalias() +=
                    L0 * f_Tv_D_Tv * dNdx.transpose() * dNdx * w / rho_LR;
                // temperature equation, pressure part
                local_K.template block<temperature_size, pressure_size>(temperature_index,
                                                       pressure_index)
                .noalias() +=
                    L0 * D_pv * dNdx.transpose() * dNdx * w / rho_LR;
            }
        }
    }

    if (_process_data.apply_mass_lumping)
    {
        auto Mpp = local_M.template block<pressure_size, pressure_size>(
            pressure_index, pressure_index);
        Mpp = Mpp.colwise().sum().eval().asDiagonal();
    }


}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
std::vector<double> const&
ThermoRichardsFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    getIntPtDarcyVelocity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    cache.clear();
    auto cache_matrix = MathLib::createZeroedMatrix<
        Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, GlobalDim, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        cache_matrix.col(ip).noalias() = _ip_data[ip].v_darcy;
    }

    return cache;
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
std::vector<double> ThermoRichardsFlowLocalAssembler<
    ShapeFunction, IntegrationMethod, GlobalDim>::getSaturation() const
{
    std::vector<double> result;
    getIntPtSaturation(0, {}, {}, result);
    return result;
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
std::vector<double> const&
ThermoRichardsFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    getIntPtSaturation(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(
        _ip_data, &IpData::saturation, cache);
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
std::vector<double> ThermoRichardsFlowLocalAssembler<
    ShapeFunction, IntegrationMethod, GlobalDim>::getPorosity() const
{
    std::vector<double> result;
    getIntPtPorosity(0, {}, {}, result);
    return result;
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
std::vector<double> const&
ThermoRichardsFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    getIntPtPorosity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(_ip_data,
                                                     &IpData::porosity, cache);
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
std::vector<double> const&
ThermoRichardsFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    getIntPtDryDensitySolid(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(
        _ip_data, &IpData::dry_density_solid, cache);
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void ThermoRichardsFlowLocalAssembler<ShapeFunction, IntegrationMethod,
                                      GlobalDim>::
    computeSecondaryVariableConcrete(double const t, double const dt,
                                     Eigen::VectorXd const& local_x,
                                     Eigen::VectorXd const& local_x_dot)
{
    auto const T =
        local_x.template segment<temperature_size>(temperature_index);

    auto const T_dot =
        local_x_dot.template segment<temperature_size>(temperature_index);

    auto const p_L = local_x.template segment<pressure_size>(pressure_index);

    auto p_L_dot = local_x_dot.template segment<pressure_size>(pressure_index);

    auto const& medium = *_process_data.media_map->getMedium(_element.getID());
    auto const& liquid_phase = medium.phase("AqueousLiquid");
    auto const& solid_phase = medium.phase("Solid");
    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    double saturation_avg = 0;
    double porosity_avg = 0;

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& N = _ip_data[ip].N;

        double T_ip;
        NumLib::shapeFunctionInterpolate(T, N, T_ip);
        double T_dot_ip;
        NumLib::shapeFunctionInterpolate(T_dot, N, T_dot_ip);

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N, p_cap_ip);

        double p_cap_dot_ip;
        NumLib::shapeFunctionInterpolate(-p_L_dot, N, p_cap_dot_ip);

        variables[static_cast<int>(MPL::Variable::capillary_pressure)] =
            p_cap_ip;
        variables[static_cast<int>(MPL::Variable::phase_pressure)] = -p_cap_ip;

        variables[static_cast<int>(MPL::Variable::temperature)] = T_ip;

        auto& S_L = _ip_data[ip].saturation;
        auto const S_L_prev = _ip_data[ip].saturation_prev;
        S_L = medium[MPL::PropertyType::saturation]
                  .template value<double>(variables, x_position, t, dt);
        variables[static_cast<int>(MPL::Variable::liquid_saturation)] = S_L;
        variables_prev[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L_prev;

        auto chi_S_L = S_L;
        auto chi_S_L_prev = S_L_prev;
        if (medium.hasProperty(MPL::PropertyType::bishops_effective_stress))
        {
            auto const chi = [&medium, x_position, t, dt](double const S_L) {
                MPL::VariableArray variables;
                variables[static_cast<int>(MPL::Variable::liquid_saturation)] =
                    S_L;
                return medium[MPL::PropertyType::bishops_effective_stress]
                    .template value<double>(variables, x_position, t, dt);
            };
            chi_S_L = chi(S_L);
            chi_S_L_prev = chi(S_L_prev);
        }
        variables[static_cast<int>(MPL::Variable::effective_pore_pressure)] =
            -chi_S_L * p_cap_ip;
        variables_prev[static_cast<int>(
            MPL::Variable::effective_pore_pressure)] =
            -chi_S_L_prev * (p_cap_ip - p_cap_dot_ip * dt);

        auto const alpha =
            medium[MPL::PropertyType::biot_coefficient]
                .template value<double>(variables, x_position, t, dt);

        auto& solid_elasticity = *_process_data.simplified_elasticity;
        auto const beta_S =
            solid_elasticity.bulkCompressibilityFromYoungsModulus(
                solid_phase, variables, x_position, t, dt);
        auto const beta_SR = (1 - alpha) * beta_S;
        variables[static_cast<int>(MPL::Variable::grain_compressibility)] =
            beta_SR;

        auto& phi = _ip_data[ip].porosity;
        {  // Porosity update
            variables_prev[static_cast<int>(MPL::Variable::porosity)] =
                _ip_data[ip].porosity_prev;
            phi = medium[MPL::PropertyType::porosity]
                      .template value<double>(variables, variables_prev,
                                              x_position, t, dt);
            variables[static_cast<int>(MPL::Variable::porosity)] = phi;
        }

        auto const mu =
            liquid_phase[MPL::PropertyType::viscosity]
                .template value<double>(variables, x_position, t, dt);
        auto const rho_LR =
            liquid_phase[MPL::PropertyType::density]
                .template value<double>(variables, x_position, t, dt);

        auto const K_intrinsic = MPL::formEigenTensor<GlobalDim>(
            medium[MPL::PropertyType::permeability]
                .value(variables, x_position, t, dt));

        double const k_rel =
            medium[MPL::PropertyType::relative_permeability]
                .template value<double>(variables, x_position, t, dt);

        GlobalDimMatrixType const K_over_mu = k_rel * K_intrinsic / mu;

        auto const rho_SR =
            solid_phase[MPL::PropertyType::density]
                .template value<double>(variables, x_position, t, dt);
        _ip_data[ip].dry_density_solid = (1 - phi) * rho_SR;

        auto const& b = _process_data.specific_body_force;

        // Compute the velocity
        auto const& dNdx = _ip_data[ip].dNdx;
        _ip_data[ip].v_darcy.noalias() =
            -K_over_mu * dNdx * p_L + rho_LR * K_over_mu * b;

        saturation_avg += S_L;
        porosity_avg += phi;
    }
    saturation_avg /= n_integration_points;
    porosity_avg /= n_integration_points;

    (*_process_data.element_saturation)[_element.getID()] = saturation_avg;
    (*_process_data.element_porosity)[_element.getID()] = porosity_avg;
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
unsigned ThermoRichardsFlowLocalAssembler<
    ShapeFunction, IntegrationMethod, GlobalDim>::getNumberOfIntegrationPoints()
    const
{
    return _integration_method.getNumberOfPoints();
}

}  // namespace ThermoRichardsFlow
}  // namespace ProcessLib
