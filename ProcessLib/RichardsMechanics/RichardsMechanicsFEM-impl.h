/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on November 29, 2017, 2:03 PM
 */

#pragma once

#include <Eigen/LU>
#include <cassert>

#include "ComputeMicroPorosity.h"
#include "IntegrationPointData.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/EigenBlockMatrixView.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "ProcessLib/Utils/TransposeInPlace.h"

namespace ProcessLib
{
namespace RichardsMechanics
{
template <int DisplacementDim, typename IPData>
void updateSwellingStressAndVolumetricStrain(
    IPData& ip_data, MaterialPropertyLib::Medium const& medium,
    MaterialPropertyLib::Phase const& solid_phase,
    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> const& C_el,
    double const rho_LR, double const mu,
    std::optional<MicroPorosityParameters> micro_porosity_parameters,
    double const alpha, double const phi, double const p_cap_ip,
    MPL::VariableArray& variables, MPL::VariableArray& variables_prev,
    ParameterLib::SpatialPosition const& x_position, double const t,
    double const dt)
{
    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    auto& sigma_sw = ip_data.sigma_sw;
    auto const& sigma_sw_prev = ip_data.sigma_sw_prev;
    if (!medium.hasProperty(MPL::PropertyType::saturation_micro))
    {
        // If there is swelling, compute it. Update volumetric strain rate,
        // s.t. it corresponds to the mechanical part only.
        sigma_sw = sigma_sw_prev;
        if (solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate))
        {
            auto const sigma_sw_dot =
                MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
                    MPL::formEigenTensor<3>(
                        solid_phase[MPL::PropertyType::swelling_stress_rate]
                            .value(variables, variables_prev, x_position, t,
                                   dt)));
            sigma_sw += sigma_sw_dot * dt;

            // !!! Misusing volumetric strain for mechanical volumetric
            // strain just to update the transport porosity !!!
            variables.volumetric_strain +=
                identity2.transpose() * C_el.inverse() * sigma_sw;
            variables_prev.volumetric_strain +=
                identity2.transpose() * C_el.inverse() * sigma_sw_prev;
        }
    }

    // TODO (naumov) saturation_micro must be always defined together with
    // the micro_porosity_parameters.
    if (medium.hasProperty(MPL::PropertyType::saturation_micro))
    {
        double const phi_M_prev = ip_data.transport_porosity_prev;
        double const phi_prev = ip_data.porosity_prev;
        double const phi_m_prev = phi_prev - phi_M_prev;
        double const p_L_m_prev = ip_data.liquid_pressure_m_prev;

        auto const S_L_m_prev = ip_data.saturation_m_prev;

        auto const [delta_phi_m, delta_e_sw, delta_p_L_m, delta_sigma_sw] =
            computeMicroPorosity<DisplacementDim>(
                identity2.transpose() * C_el.inverse(), rho_LR, mu,
                *micro_porosity_parameters, alpha, phi, -p_cap_ip, p_L_m_prev,
                variables_prev, S_L_m_prev, phi_m_prev, x_position, t, dt,
                medium.property(MPL::PropertyType::saturation_micro),
                solid_phase.property(MPL::PropertyType::swelling_stress_rate));

        auto& phi_M = ip_data.transport_porosity;
        phi_M = phi - (phi_m_prev + delta_phi_m);
        variables_prev.transport_porosity = phi_M_prev;
        variables.transport_porosity = phi_M;

        auto& p_L_m = ip_data.liquid_pressure_m;
        p_L_m = p_L_m_prev + delta_p_L_m;
        {  // Update micro saturation.
            MPL::VariableArray variables_prev;
            variables_prev.capillary_pressure = -p_L_m_prev;
            MPL::VariableArray variables;
            variables.capillary_pressure = -p_L_m;

            ip_data.saturation_m =
                medium.property(MPL::PropertyType::saturation_micro)
                    .template value<double>(variables, x_position, t, dt);
        }
        sigma_sw = sigma_sw_prev + delta_sigma_sw;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                ShapeFunctionPressure, DisplacementDim>::
    RichardsMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        RichardsMechanicsProcessData<DisplacementDim>& process_data)
    : _process_data(process_data),
      _integration_method(integration_method),
      _element(e),
      _is_axially_symmetric(is_axially_symmetric)
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);
    _secondary_data.N_u.resize(n_integration_points);

    auto const shape_matrices_u =
        NumLib::initShapeMatrices<ShapeFunctionDisplacement,
                                  ShapeMatricesTypeDisplacement,
                                  DisplacementDim>(e, is_axially_symmetric,
                                                   _integration_method);

    auto const shape_matrices_p =
        NumLib::initShapeMatrices<ShapeFunctionPressure,
                                  ShapeMatricesTypePressure, DisplacementDim>(
            e, is_axially_symmetric, _integration_method);

    auto const& solid_material =
        MaterialLib::Solids::selectSolidConstitutiveRelation(
            _process_data.solid_materials, _process_data.material_ids,
            e.getID());

    auto const& medium = _process_data.media_map.getMedium(_element.getID());

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        _ip_data.emplace_back(solid_material);
        auto& ip_data = _ip_data[ip];
        auto const& sm_u = shape_matrices_u[ip];
        _ip_data[ip].integration_weight =
            _integration_method.getWeightedPoint(ip).getWeight() *
            sm_u.integralMeasure * sm_u.detJ;

        ip_data.N_u = sm_u.N;
        ip_data.dNdx_u = sm_u.dNdx;

        ip_data.N_p = shape_matrices_p[ip].N;
        ip_data.dNdx_p = shape_matrices_p[ip].dNdx;

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

        _secondary_data.N_u[ip] = shape_matrices_u[ip].N;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::size_t RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::setIPDataInitialConditions(std::string_view name,
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

    if (name == "sigma")
    {
        if (_process_data.initial_stress != nullptr)
        {
            OGS_FATAL(
                "Setting initial conditions for stress from integration "
                "point data and from a parameter '{:s}' is not possible "
                "simultaneously.",
                _process_data.initial_stress->name);
        }
        return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
            values, _ip_data, &IpData::sigma_eff);
    }

    if (name == "saturation")
    {
        return ProcessLib::setIntegrationPointScalarData(values, _ip_data,
                                                         &IpData::saturation);
    }
    if (name == "porosity")
    {
        return ProcessLib::setIntegrationPointScalarData(values, _ip_data,
                                                         &IpData::porosity);
    }
    if (name == "transport_porosity")
    {
        return ProcessLib::setIntegrationPointScalarData(
            values, _ip_data, &IpData::transport_porosity);
    }
    if (name == "swelling_stress")
    {
        return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
            values, _ip_data, &IpData::sigma_sw);
    }
    if (name == "epsilon")
    {
        return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
            values, _ip_data, &IpData::eps);
    }
    if (name.starts_with("material_state_variable_"))
    {
        name.remove_prefix(24);

        // Using first ip data for solid material. TODO (naumov) move solid
        // material into element, store only material state in IPs.
        auto const& internal_variables =
            _ip_data[0].solid_material.getInternalVariables();
        if (auto const iv = std::find_if(
                begin(internal_variables), end(internal_variables),
                [&name](auto const& iv) { return iv.name == name; });
            iv != end(internal_variables))
        {
            DBUG("Setting material state variable '{:s}'", name);
            return ProcessLib::setIntegrationPointDataMaterialStateVariables(
                values, _ip_data, &IpData::material_state_variables,
                iv->reference);
        }

        ERR("Could not find variable {:s} in solid material model's internal "
            "variables.",
            name);
    }
    return 0;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, DisplacementDim>::
    setInitialConditionsConcrete(Eigen::VectorXd const local_x,
                                 double const t,
                                 int const /*process_id*/)
{
    assert(local_x.size() == pressure_size + displacement_size);

    auto const p_L = local_x.template segment<pressure_size>(pressure_index);

    constexpr double dt = std::numeric_limits<double>::quiet_NaN();
    auto const& medium = _process_data.media_map.getMedium(_element.getID());
    MPL::VariableArray variables;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    auto const& solid_phase = medium->phase("Solid");

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto const& N_p = _ip_data[ip].N_p;

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N_p, p_cap_ip);

        variables.capillary_pressure = p_cap_ip;
        variables.liquid_phase_pressure = -p_cap_ip;
        // setting pG to 1 atm
        // TODO : rewrite equations s.t. p_L = pG-p_cap
        variables.gas_phase_pressure = 1.0e5;
        _ip_data[ip].liquid_pressure_m_prev = -p_cap_ip;
        _ip_data[ip].liquid_pressure_m = -p_cap_ip;

        auto const temperature =
            medium->property(MPL::PropertyType::reference_temperature)
                .template value<double>(variables, x_position, t, dt);
        variables.temperature = temperature;

        _ip_data[ip].saturation_prev =
            medium->property(MPL::PropertyType::saturation)
                .template value<double>(variables, x_position, t, dt);

        if (medium->hasProperty(MPL::PropertyType::saturation_micro))
        {
            MPL::VariableArray vars;
            vars.capillary_pressure = p_cap_ip;
            double const S_L_m =
                medium->property(MPL::PropertyType::saturation_micro)
                    .template value<double>(vars, x_position, t, dt);
            _ip_data[ip].saturation_m_prev = S_L_m;
        }

        // Set eps_m_prev from potentially non-zero eps and sigma_sw from
        // restart.
        auto const C_el = _ip_data[ip].computeElasticTangentStiffness(
            t, x_position, dt, temperature);
        auto& eps = _ip_data[ip].eps;
        auto& sigma_sw = _ip_data[ip].sigma_sw;

        _ip_data[ip].eps_m_prev.noalias() =
            solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate)
                ? eps + C_el.inverse() * sigma_sw
                : eps;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::assemble(double const t, double const dt,
                               std::vector<double> const& local_x,
                               std::vector<double> const& local_x_prev,
                               std::vector<double>& local_M_data,
                               std::vector<double>& local_K_data,
                               std::vector<double>& local_rhs_data)
{
    assert(local_x.size() == pressure_size + displacement_size);

    auto p_L =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            pressure_size> const>(local_x.data() + pressure_index,
                                  pressure_size);

    auto u =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x.data() + displacement_index,
                                      displacement_size);

    auto p_L_prev =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            pressure_size> const>(local_x_prev.data() + pressure_index,
                                  pressure_size);

    auto u_prev =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x_prev.data() + displacement_index,
                                      displacement_size);

    auto K = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size + pressure_size,
            displacement_size + pressure_size>>(
        local_K_data, displacement_size + pressure_size,
        displacement_size + pressure_size);

    auto M = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size + pressure_size,
            displacement_size + pressure_size>>(
        local_M_data, displacement_size + pressure_size,
        displacement_size + pressure_size);

    auto rhs = MathLib::createZeroedVector<
        typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size + pressure_size>>(
        local_rhs_data, displacement_size + pressure_size);

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    auto const& medium = _process_data.media_map.getMedium(_element.getID());
    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& solid_phase = medium->phase("Solid");
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

        auto const& N_u = _ip_data[ip].N_u;
        auto const& dNdx_u = _ip_data[ip].dNdx_u;

        auto const& N_p = _ip_data[ip].N_p;
        auto const& dNdx_p = _ip_data[ip].dNdx_p;

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                _element, N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        auto& eps = _ip_data[ip].eps;
        auto& eps_m = _ip_data[ip].eps_m;
        eps.noalias() = B * u;

        auto& S_L = _ip_data[ip].saturation;
        auto const S_L_prev = _ip_data[ip].saturation_prev;

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N_p, p_cap_ip);

        double p_cap_prev_ip;
        NumLib::shapeFunctionInterpolate(-p_L_prev, N_p, p_cap_prev_ip);

        variables.capillary_pressure = p_cap_ip;
        variables.liquid_phase_pressure = -p_cap_ip;
        // setting pG to 1 atm
        // TODO : rewrite equations s.t. p_L = pG-p_cap
        variables.gas_phase_pressure = 1.0e5;

        auto const temperature =
            medium->property(MPL::PropertyType::reference_temperature)
                .template value<double>(variables, x_position, t, dt);
        variables.temperature = temperature;

        auto const alpha =
            medium->property(MPL::PropertyType::biot_coefficient)
                .template value<double>(variables, x_position, t, dt);
        auto const C_el = _ip_data[ip].computeElasticTangentStiffness(
            t, x_position, dt, temperature);

        auto const beta_SR =
            (1 - alpha) /
            _ip_data[ip].solid_material.getBulkModulus(t, x_position, &C_el);
        variables.grain_compressibility = beta_SR;

        auto const rho_LR =
            liquid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);
        variables.density = rho_LR;

        auto const& b = _process_data.specific_body_force;

        S_L = medium->property(MPL::PropertyType::saturation)
                  .template value<double>(variables, x_position, t, dt);
        variables.liquid_saturation = S_L;
        variables_prev.liquid_saturation = S_L_prev;

        // tangent derivative for Jacobian
        double const dS_L_dp_cap =
            medium->property(MPL::PropertyType::saturation)
                .template dValue<double>(variables,
                                         MPL::Variable::capillary_pressure,
                                         x_position, t, dt);
        // secant derivative from time discretization for storage
        // use tangent, if secant is not available
        double const DeltaS_L_Deltap_cap =
            (p_cap_ip == p_cap_prev_ip)
                ? dS_L_dp_cap
                : (S_L - S_L_prev) / (p_cap_ip - p_cap_prev_ip);

        auto const chi = [medium, x_position, t, dt](double const S_L)
        {
            MPL::VariableArray vs;
            vs.liquid_saturation = S_L;
            return medium->property(MPL::PropertyType::bishops_effective_stress)
                .template value<double>(vs, x_position, t, dt);
        };
        double const chi_S_L = chi(S_L);
        double const chi_S_L_prev = chi(S_L_prev);

        double const p_FR = -chi_S_L * p_cap_ip;
        variables.effective_pore_pressure = p_FR;
        variables_prev.effective_pore_pressure = -chi_S_L_prev * p_cap_prev_ip;

        // Set volumetric strain rate for the general case without swelling.
        variables.volumetric_strain = Invariants::trace(_ip_data[ip].eps);
        variables_prev.volumetric_strain = Invariants::trace(B * u_prev);

        auto& phi = _ip_data[ip].porosity;
        {  // Porosity update
            variables_prev.porosity = _ip_data[ip].porosity_prev;
            phi = medium->property(MPL::PropertyType::porosity)
                      .template value<double>(variables, variables_prev,
                                              x_position, t, dt);
            variables.porosity = phi;
        }

        if (alpha < phi)
        {
            OGS_FATAL(
                "RichardsMechanics: Biot-coefficient {} is smaller than "
                "porosity {} in element/integration point {}/{}.",
                alpha, phi, _element.getID(), ip);
        }

        // Swelling and possibly volumetric strain rate update.
        {
            auto& sigma_sw = _ip_data[ip].sigma_sw;
            auto const& sigma_sw_prev = _ip_data[ip].sigma_sw_prev;

            // If there is swelling, compute it. Update volumetric strain rate,
            // s.t. it corresponds to the mechanical part only.
            sigma_sw = sigma_sw_prev;
            if (solid_phase.hasProperty(
                    MPL::PropertyType::swelling_stress_rate))
            {
                auto const sigma_sw_dot =
                    MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
                        MPL::formEigenTensor<3>(
                            solid_phase[MPL::PropertyType::swelling_stress_rate]
                                .value(variables, variables_prev, x_position, t,
                                       dt)));
                sigma_sw += sigma_sw_dot * dt;

                // !!! Misusing volumetric strain for mechanical volumetric
                // strain just to update the transport porosity !!!
                variables.volumetric_strain +=
                    identity2.transpose() * C_el.inverse() * sigma_sw;
                variables_prev.volumetric_strain +=
                    identity2.transpose() * C_el.inverse() * sigma_sw_prev;
            }

            if (medium->hasProperty(MPL::PropertyType::transport_porosity))
            {
                variables_prev.transport_porosity =
                    _ip_data[ip].transport_porosity_prev;

                _ip_data[ip].transport_porosity =
                    medium->property(MPL::PropertyType::transport_porosity)
                        .template value<double>(variables, variables_prev,
                                                x_position, t, dt);
                variables.transport_porosity = _ip_data[ip].transport_porosity;
            }
            else
            {
                variables.transport_porosity = phi;
            }
        }

        double const k_rel =
            medium->property(MPL::PropertyType::relative_permeability)
                .template value<double>(variables, x_position, t, dt);
        auto const mu =
            liquid_phase.property(MPL::PropertyType::viscosity)
                .template value<double>(variables, x_position, t, dt);

        auto const& sigma_sw = _ip_data[ip].sigma_sw;
        auto const& sigma_eff = _ip_data[ip].sigma_eff;

        // Set mechanical variables for the intrinsic permeability model
        // For stress dependent permeability.
        {
            auto const sigma_total =
                (sigma_eff - alpha * p_FR * identity2).eval();

            // For stress dependent permeability.
            variables.total_stress.emplace<SymmetricTensor>(
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                    sigma_total));
        }

        variables.equivalent_plastic_strain =
            _ip_data[ip].material_state_variables->getEquivalentPlasticStrain();

        auto const K_intrinsic = MPL::formEigenTensor<DisplacementDim>(
            medium->property(MPL::PropertyType::permeability)
                .value(variables, x_position, t, dt));

        GlobalDimMatrixType const rho_K_over_mu =
            K_intrinsic * rho_LR * k_rel / mu;

        //
        // displacement equation, displacement part
        //
        eps_m.noalias() =
            solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate)
                ? eps + C_el.inverse() * sigma_sw
                : eps;
        variables.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m);

        _ip_data[ip].updateConstitutiveRelation(variables, t, x_position, dt,
                                                temperature);

        // p_SR
        variables.solid_grain_pressure =
            p_FR - sigma_eff.dot(identity2) / (3 * (1 - phi));
        auto const rho_SR =
            solid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);

        //
        // displacement equation, displacement part
        //
        double const rho = rho_SR * (1 - phi) + S_L * phi * rho_LR;
        rhs.template segment<displacement_size>(displacement_index).noalias() -=
            (B.transpose() * sigma_eff - N_u_op(N_u).transpose() * rho * b) * w;

        //
        // pressure equation, pressure part.
        //
        auto const beta_LR =
            1 / rho_LR *
            liquid_phase.property(MPL::PropertyType::density)
                .template dValue<double>(variables,
                                         MPL::Variable::liquid_phase_pressure,
                                         x_position, t, dt);

        double const a0 = S_L * (alpha - phi) * beta_SR;
        // Volumetric average specific storage of the solid and fluid phases.
        double const specific_storage =
            DeltaS_L_Deltap_cap * (p_cap_ip * a0 - phi) +
            S_L * (phi * beta_LR + a0);
        M.template block<pressure_size, pressure_size>(pressure_index,
                                                       pressure_index)
            .noalias() += N_p.transpose() * rho_LR * specific_storage * N_p * w;

        K.template block<pressure_size, pressure_size>(pressure_index,
                                                       pressure_index)
            .noalias() += dNdx_p.transpose() * rho_K_over_mu * dNdx_p * w;

        rhs.template segment<pressure_size>(pressure_index).noalias() +=
            dNdx_p.transpose() * rho_LR * rho_K_over_mu * b * w;

        //
        // displacement equation, pressure part
        //
        K.template block<displacement_size, pressure_size>(displacement_index,
                                                           pressure_index)
            .noalias() -= B.transpose() * alpha * chi_S_L * identity2 * N_p * w;

        //
        // pressure equation, displacement part.
        //
        M.template block<pressure_size, displacement_size>(pressure_index,
                                                           displacement_index)
            .noalias() += N_p.transpose() * S_L * rho_LR * alpha *
                          identity2.transpose() * B * w;
    }

    if (_process_data.apply_mass_lumping)
    {
        auto Mpp = M.template block<pressure_size, pressure_size>(
            pressure_index, pressure_index);
        Mpp = Mpp.colwise().sum().eval().asDiagonal();
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, DisplacementDim>::
    assembleWithJacobian(double const t, double const dt,
                         std::vector<double> const& local_x,
                         std::vector<double> const& local_x_prev,
                         std::vector<double>& /*local_M_data*/,
                         std::vector<double>& /*local_K_data*/,
                         std::vector<double>& local_rhs_data,
                         std::vector<double>& local_Jac_data)
{
    assert(local_x.size() == pressure_size + displacement_size);

    auto p_L =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            pressure_size> const>(local_x.data() + pressure_index,
                                  pressure_size);

    auto u =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x.data() + displacement_index,
                                      displacement_size);

    auto p_L_prev =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            pressure_size> const>(local_x_prev.data() + pressure_index,
                                  pressure_size);
    auto u_prev =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x_prev.data() + displacement_index,
                                      displacement_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size + pressure_size,
            displacement_size + pressure_size>>(
        local_Jac_data, displacement_size + pressure_size,
        displacement_size + pressure_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size + pressure_size>>(
        local_rhs_data, displacement_size + pressure_size);

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    typename ShapeMatricesTypePressure::NodalMatrixType laplace_p =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage_p_a_p =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage_p_a_S_Jpp =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage_p_a_S =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        displacement_size, pressure_size>
        Kup = ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size, pressure_size>::Zero(displacement_size,
                                                    pressure_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        pressure_size, displacement_size>
        Kpu = ShapeMatricesTypeDisplacement::template MatrixType<
            pressure_size, displacement_size>::Zero(pressure_size,
                                                    displacement_size);

    auto const& medium = _process_data.media_map.getMedium(_element.getID());
    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& solid_phase = medium->phase("Solid");
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

        auto const& N_u = _ip_data[ip].N_u;
        auto const& dNdx_u = _ip_data[ip].dNdx_u;

        auto const& N_p = _ip_data[ip].N_p;
        auto const& dNdx_p = _ip_data[ip].dNdx_p;

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                _element, N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N_p, p_cap_ip);

        double p_cap_prev_ip;
        NumLib::shapeFunctionInterpolate(-p_L_prev, N_p, p_cap_prev_ip);

        variables.capillary_pressure = p_cap_ip;
        variables.liquid_phase_pressure = -p_cap_ip;
        // setting pG to 1 atm
        // TODO : rewrite equations s.t. p_L = pG-p_cap
        variables.gas_phase_pressure = 1.0e5;

        auto const temperature =
            medium->property(MPL::PropertyType::reference_temperature)
                .template value<double>(variables, x_position, t, dt);
        variables.temperature = temperature;

        auto& eps = _ip_data[ip].eps;
        auto& eps_m = _ip_data[ip].eps_m;
        eps.noalias() = B * u;
        auto const& sigma_eff = _ip_data[ip].sigma_eff;
        auto& S_L = _ip_data[ip].saturation;
        auto const S_L_prev = _ip_data[ip].saturation_prev;
        auto const alpha =
            medium->property(MPL::PropertyType::biot_coefficient)
                .template value<double>(variables, x_position, t, dt);

        auto const C_el = _ip_data[ip].computeElasticTangentStiffness(
            t, x_position, dt, temperature);

        auto const beta_SR =
            (1 - alpha) /
            _ip_data[ip].solid_material.getBulkModulus(t, x_position, &C_el);
        variables.grain_compressibility = beta_SR;

        auto const rho_LR =
            liquid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);
        variables.density = rho_LR;

        auto const& b = _process_data.specific_body_force;

        S_L = medium->property(MPL::PropertyType::saturation)
                  .template value<double>(variables, x_position, t, dt);
        variables.liquid_saturation = S_L;
        variables_prev.liquid_saturation = S_L_prev;

        // tangent derivative for Jacobian
        double const dS_L_dp_cap =
            medium->property(MPL::PropertyType::saturation)
                .template dValue<double>(variables,
                                         MPL::Variable::capillary_pressure,
                                         x_position, t, dt);
        // secant derivative from time discretization for storage
        // use tangent, if secant is not available
        double const DeltaS_L_Deltap_cap =
            (p_cap_ip == p_cap_prev_ip)
                ? dS_L_dp_cap
                : (S_L - S_L_prev) / (p_cap_ip - p_cap_prev_ip);

        auto const chi = [medium, x_position, t, dt](double const S_L)
        {
            MPL::VariableArray vs;
            vs.liquid_saturation = S_L;
            return medium->property(MPL::PropertyType::bishops_effective_stress)
                .template value<double>(vs, x_position, t, dt);
        };
        double const chi_S_L = chi(S_L);
        double const chi_S_L_prev = chi(S_L_prev);

        double const p_FR = -chi_S_L * p_cap_ip;
        variables.effective_pore_pressure = p_FR;
        variables_prev.effective_pore_pressure = -chi_S_L_prev * p_cap_prev_ip;

        // Set volumetric strain rate for the general case without swelling.
        variables.volumetric_strain = Invariants::trace(_ip_data[ip].eps);
        variables_prev.volumetric_strain = Invariants::trace(B * u_prev);

        auto& phi = _ip_data[ip].porosity;
        {  // Porosity update

            variables_prev.porosity = _ip_data[ip].porosity_prev;
            phi = medium->property(MPL::PropertyType::porosity)
                      .template value<double>(variables, variables_prev,
                                              x_position, t, dt);
            variables.porosity = phi;
        }

        if (alpha < phi)
        {
            OGS_FATAL(
                "RichardsMechanics: Biot-coefficient {} is smaller than "
                "porosity {} in element/integration point {}/{}.",
                alpha, phi, _element.getID(), ip);
        }

        auto const mu =
            liquid_phase.property(MPL::PropertyType::viscosity)
                .template value<double>(variables, x_position, t, dt);

        // Swelling and possibly volumetric strain rate update.
        updateSwellingStressAndVolumetricStrain<DisplacementDim>(
            _ip_data[ip], *medium, solid_phase, C_el, rho_LR, mu,
            _process_data.micro_porosity_parameters, alpha, phi, p_cap_ip,
            variables, variables_prev, x_position, t, dt);

        if (medium->hasProperty(MPL::PropertyType::transport_porosity))
        {
            if (!medium->hasProperty(MPL::PropertyType::saturation_micro))
            {
                variables_prev.transport_porosity =
                    _ip_data[ip].transport_porosity_prev;

                _ip_data[ip].transport_porosity =
                    medium->property(MPL::PropertyType::transport_porosity)
                        .template value<double>(variables, variables_prev,
                                                x_position, t, dt);
                variables.transport_porosity = _ip_data[ip].transport_porosity;
            }
        }
        else
        {
            variables.transport_porosity = phi;
        }

        double const k_rel =
            medium->property(MPL::PropertyType::relative_permeability)
                .template value<double>(variables, x_position, t, dt);

        // Set mechanical variables for the intrinsic permeability model
        // For stress dependent permeability.
        {
            auto const sigma_total =
                (_ip_data[ip].sigma_eff + alpha * p_FR * identity2).eval();
            // For stress dependent permeability.
            variables.total_stress.emplace<SymmetricTensor>(
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                    sigma_total));
        }

        variables.equivalent_plastic_strain =
            _ip_data[ip].material_state_variables->getEquivalentPlasticStrain();

        auto const K_intrinsic = MPL::formEigenTensor<DisplacementDim>(
            medium->property(MPL::PropertyType::permeability)
                .value(variables, x_position, t, dt));

        GlobalDimMatrixType const rho_Ki_over_mu = K_intrinsic * rho_LR / mu;

        //
        // displacement equation, displacement part
        //
        auto& sigma_sw = _ip_data[ip].sigma_sw;

        eps_m.noalias() =
            solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate)
                ? eps + C_el.inverse() * sigma_sw
                : eps;
        variables.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m);

        auto C = _ip_data[ip].updateConstitutiveRelation(
            variables, t, x_position, dt, temperature);

        local_Jac
            .template block<displacement_size, displacement_size>(
                displacement_index, displacement_index)
            .noalias() += B.transpose() * C * B * w;

        // p_SR
        variables.solid_grain_pressure =
            p_FR - sigma_eff.dot(identity2) / (3 * (1 - phi));
        auto const rho_SR =
            solid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);

        double const rho = rho_SR * (1 - phi) + S_L * phi * rho_LR;
        local_rhs.template segment<displacement_size>(displacement_index)
            .noalias() -=
            (B.transpose() * sigma_eff - N_u_op(N_u).transpose() * rho * b) * w;

        //
        // displacement equation, pressure part
        //
        Kup.noalias() += B.transpose() * alpha * chi_S_L * identity2 * N_p * w;

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
                          identity2 * N_p * w;

        local_Jac
            .template block<displacement_size, pressure_size>(
                displacement_index, pressure_index)
            .noalias() +=
            N_u_op(N_u).transpose() * phi * rho_LR * dS_L_dp_cap * b * N_p * w;

        // For the swelling stress with double structure model the corresponding
        // Jacobian u-p entry would be required, but it does not improve
        // convergence and sometimes worsens it:
        // if (medium->hasProperty(MPL::PropertyType::saturation_micro))
        // {
        //     -B.transpose() *
        //         dsigma_sw_dS_L_m* dS_L_m_dp_cap_m*(p_L_m - p_L_m_prev) /
        //         (p_cap_ip - p_cap_prev_ip) * N_p* w;
        // }
        if (!medium->hasProperty(MPL::PropertyType::saturation_micro) &&
            solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate))
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
                B.transpose() * dsigma_sw_dS_L * dS_L_dp_cap * N_p * w;
        }
        //
        // pressure equation, displacement part.
        //
        if (_process_data.explicit_hm_coupling_in_unsaturated_zone)
        {
            Kpu.noalias() += N_p.transpose() * chi_S_L_prev * rho_LR * alpha *
                             identity2.transpose() * B * w;
        }
        else
        {
            Kpu.noalias() += N_p.transpose() * S_L * rho_LR * alpha *
                             identity2.transpose() * B * w;
        }

        //
        // pressure equation, pressure part.
        //
        laplace_p.noalias() +=
            dNdx_p.transpose() * k_rel * rho_Ki_over_mu * dNdx_p * w;

        auto const beta_LR =
            1 / rho_LR *
            liquid_phase.property(MPL::PropertyType::density)
                .template dValue<double>(variables,
                                         MPL::Variable::liquid_phase_pressure,
                                         x_position, t, dt);

        double const a0 = (alpha - phi) * beta_SR;
        double const specific_storage_a_p = S_L * (phi * beta_LR + S_L * a0);
        double const specific_storage_a_S = phi - p_cap_ip * S_L * a0;

        double const dspecific_storage_a_p_dp_cap =
            dS_L_dp_cap * (phi * beta_LR + 2 * S_L * a0);
        double const dspecific_storage_a_S_dp_cap =
            -a0 * (S_L + p_cap_ip * dS_L_dp_cap);

        storage_p_a_p.noalias() +=
            N_p.transpose() * rho_LR * specific_storage_a_p * N_p * w;

        storage_p_a_S.noalias() -= N_p.transpose() * rho_LR *
                                   specific_storage_a_S * DeltaS_L_Deltap_cap *
                                   N_p * w;

        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() += N_p.transpose() * (p_cap_ip - p_cap_prev_ip) / dt *
                          rho_LR * dspecific_storage_a_p_dp_cap * N_p * w;

        storage_p_a_S_Jpp.noalias() -=
            N_p.transpose() * rho_LR *
            ((S_L - S_L_prev) * dspecific_storage_a_S_dp_cap +
             specific_storage_a_S * dS_L_dp_cap) /
            dt * N_p * w;

        if (!_process_data.explicit_hm_coupling_in_unsaturated_zone)
        {
            local_Jac
                .template block<pressure_size, pressure_size>(pressure_index,
                                                              pressure_index)
                .noalias() -= N_p.transpose() * rho_LR * dS_L_dp_cap * alpha *
                              identity2.transpose() * B * (u - u_prev) / dt *
                              N_p * w;
        }

        double const dk_rel_dS_l =
            medium->property(MPL::PropertyType::relative_permeability)
                .template dValue<double>(variables,
                                         MPL::Variable::liquid_saturation,
                                         x_position, t, dt);
        typename ShapeMatricesTypeDisplacement::GlobalDimVectorType const
            grad_p_cap = -dNdx_p * p_L;
        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() += dNdx_p.transpose() * rho_Ki_over_mu * grad_p_cap *
                          dk_rel_dS_l * dS_L_dp_cap * N_p * w;

        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() += dNdx_p.transpose() * rho_LR * rho_Ki_over_mu * b *
                          dk_rel_dS_l * dS_L_dp_cap * N_p * w;

        local_rhs.template segment<pressure_size>(pressure_index).noalias() +=
            dNdx_p.transpose() * rho_LR * k_rel * rho_Ki_over_mu * b * w;

        if (medium->hasProperty(MPL::PropertyType::saturation_micro))
        {
            double const alpha_bar = _process_data.micro_porosity_parameters
                                         ->mass_exchange_coefficient;
            auto const p_L_m = _ip_data[ip].liquid_pressure_m;
            local_rhs.template segment<pressure_size>(pressure_index)
                .noalias() -=
                N_p.transpose() * alpha_bar / mu * (-p_cap_ip - p_L_m) * w;

            local_Jac
                .template block<pressure_size, pressure_size>(pressure_index,
                                                              pressure_index)
                .noalias() += N_p.transpose() * alpha_bar / mu * N_p * w;
            if (p_cap_ip != p_cap_prev_ip)
            {
                double const p_L_m_prev = _ip_data[ip].liquid_pressure_m_prev;
                local_Jac
                    .template block<pressure_size, pressure_size>(
                        pressure_index, pressure_index)
                    .noalias() += N_p.transpose() * alpha_bar / mu *
                                  (p_L_m - p_L_m_prev) /
                                  (p_cap_ip - p_cap_prev_ip) * N_p * w;
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

    // pressure equation, pressure part.
    local_Jac
        .template block<pressure_size, pressure_size>(pressure_index,
                                                      pressure_index)
        .noalias() += laplace_p + storage_p_a_p / dt + storage_p_a_S_Jpp;

    // pressure equation, displacement part.
    local_Jac
        .template block<pressure_size, displacement_size>(pressure_index,
                                                          displacement_index)
        .noalias() = Kpu / dt;

    // pressure equation
    local_rhs.template segment<pressure_size>(pressure_index).noalias() -=
        laplace_p * p_L +
        (storage_p_a_p + storage_p_a_S) * (p_L - p_L_prev) / dt +
        Kpu * (u - u_prev) / dt;

    // displacement equation
    local_rhs.template segment<displacement_size>(displacement_index)
        .noalias() += Kup * p_L;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::getSigma() const
{
    constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);

    return transposeInPlace<kelvin_vector_size>(
        [this](std::vector<double>& values)
        { return getIntPtSigma(0, {}, {}, values); });
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> const& RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    getIntPtSigma(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
        _ip_data, &IpData::sigma_eff, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::getSwellingStress() const
{
    constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);

    return transposeInPlace<kelvin_vector_size>(
        [this](std::vector<double>& values)
        { return getIntPtSwellingStress(0, {}, {}, values); });
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> const& RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    getIntPtSwellingStress(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    auto const n_integration_points = _ip_data.size();

    cache.clear();
    auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, kelvin_vector_size, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        auto const& sigma_sw = _ip_data[ip].sigma_sw;
        cache_mat.col(ip) =
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_sw);
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::getEpsilon() const
{
    constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);

    return transposeInPlace<kelvin_vector_size>(
        [this](std::vector<double>& values)
        { return getIntPtEpsilon(0, {}, {}, values); });
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> const& RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    getIntPtEpsilon(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
        _ip_data, &IpData::eps, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
int RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                    ShapeFunctionPressure,
                                    DisplacementDim>::getMaterialID() const
{
    return _process_data.material_ids == nullptr
               ? 0
               : (*_process_data.material_ids)[_element.getID()];
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    getMaterialStateVariableInternalState(
        std::function<std::span<double>(
            typename MaterialLib::Solids::MechanicsBase<DisplacementDim>::
                MaterialStateVariables&)> const& get_values_span,
        int const& n_components) const
{
    return ProcessLib::getIntegrationPointDataMaterialStateVariables(
        _ip_data, &IpData::material_state_variables, get_values_span,
        n_components);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> const& RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    getIntPtDarcyVelocity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    cache.clear();
    auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, DisplacementDim, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        cache_matrix.col(ip).noalias() = _ip_data[ip].v_darcy;
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::getSaturation() const
{
    std::vector<double> result;
    getIntPtSaturation(0, {}, {}, result);
    return result;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> const& RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    getIntPtSaturation(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(
        _ip_data, &IpData::saturation, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::getMicroSaturation() const
{
    std::vector<double> result;
    getIntPtMicroSaturation(0, {}, {}, result);
    return result;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> const& RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    getIntPtMicroSaturation(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(
        _ip_data, &IpData::saturation_m, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::getMicroPressure() const
{
    std::vector<double> result;
    getIntPtMicroPressure(0, {}, {}, result);
    return result;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> const& RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    getIntPtMicroPressure(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(
        _ip_data, &IpData::liquid_pressure_m, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::getPorosity() const
{
    std::vector<double> result;
    getIntPtPorosity(0, {}, {}, result);
    return result;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> const& RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    getIntPtPorosity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(_ip_data,
                                                     &IpData::porosity, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::getTransportPorosity() const
{
    std::vector<double> result;
    getIntPtTransportPorosity(0, {}, {}, result);
    return result;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> const& RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    getIntPtTransportPorosity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(
        _ip_data, &IpData::transport_porosity, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> const& RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    getIntPtDryDensitySolid(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(
        _ip_data, &IpData::dry_density_solid, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, DisplacementDim>::
    assembleWithJacobianForPressureEquations(
        const double /*t*/, double const /*dt*/,
        Eigen::VectorXd const& /*local_x*/,
        Eigen::VectorXd const& /*local_x_prev*/,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& /*local_b_data*/,
        std::vector<double>& /*local_Jac_data*/)
{
    OGS_FATAL("RichardsMechanics; The staggered scheme is not implemented.");
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, DisplacementDim>::
    assembleWithJacobianForDeformationEquations(
        const double /*t*/, double const /*dt*/,
        Eigen::VectorXd const& /*local_x*/,
        Eigen::VectorXd const& /*local_x_prev*/,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& /*local_b_data*/,
        std::vector<double>& /*local_Jac_data*/)
{
    OGS_FATAL("RichardsMechanics; The staggered scheme is not implemented.");
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, DisplacementDim>::
    assembleWithJacobianForStaggeredScheme(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_x_prev, int const process_id,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    // For the equations with pressure
    if (process_id == 0)
    {
        assembleWithJacobianForPressureEquations(t, dt, local_x, local_x_prev,
                                                 local_M_data, local_K_data,
                                                 local_b_data, local_Jac_data);
        return;
    }

    // For the equations with deformation
    assembleWithJacobianForDeformationEquations(t, dt, local_x, local_x_prev,
                                                local_M_data, local_K_data,
                                                local_b_data, local_Jac_data);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, DisplacementDim>::
    computeSecondaryVariableConcrete(double const t, double const dt,
                                     Eigen::VectorXd const& local_x,
                                     Eigen::VectorXd const& local_x_prev)
{
    auto p_L = local_x.template segment<pressure_size>(pressure_index);
    auto u = local_x.template segment<displacement_size>(displacement_index);

    auto p_L_prev =
        local_x_prev.template segment<pressure_size>(pressure_index);
    auto u_prev =
        local_x_prev.template segment<displacement_size>(displacement_index);

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    auto const& medium = _process_data.media_map.getMedium(_element.getID());
    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& solid_phase = medium->phase("Solid");
    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    double saturation_avg = 0;
    double porosity_avg = 0;

    using KV = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    KV sigma_avg = KV::Zero();

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& N_p = _ip_data[ip].N_p;
        auto const& N_u = _ip_data[ip].N_u;
        auto const& dNdx_u = _ip_data[ip].dNdx_u;

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                _element, N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N_p, p_cap_ip);

        double p_cap_prev_ip;
        NumLib::shapeFunctionInterpolate(-p_L_prev, N_p, p_cap_prev_ip);

        variables.capillary_pressure = p_cap_ip;
        variables.liquid_phase_pressure = -p_cap_ip;
        // setting pG to 1 atm
        // TODO : rewrite equations s.t. p_L = pG-p_cap
        variables.gas_phase_pressure = 1.0e5;

        auto const temperature =
            medium->property(MPL::PropertyType::reference_temperature)
                .template value<double>(variables, x_position, t, dt);
        variables.temperature = temperature;

        auto& eps = _ip_data[ip].eps;
        eps.noalias() = B * u;
        auto& eps_m = _ip_data[ip].eps_m;
        auto& S_L = _ip_data[ip].saturation;
        auto const S_L_prev = _ip_data[ip].saturation_prev;
        S_L = medium->property(MPL::PropertyType::saturation)
                  .template value<double>(variables, x_position, t, dt);
        variables.liquid_saturation = S_L;
        variables_prev.liquid_saturation = S_L_prev;

        auto const chi = [medium, x_position, t, dt](double const S_L)
        {
            MPL::VariableArray vs;
            vs.liquid_saturation = S_L;
            return medium->property(MPL::PropertyType::bishops_effective_stress)
                .template value<double>(vs, x_position, t, dt);
        };
        double const chi_S_L = chi(S_L);
        double const chi_S_L_prev = chi(S_L_prev);

        auto const alpha =
            medium->property(MPL::PropertyType::biot_coefficient)
                .template value<double>(variables, x_position, t, dt);

        auto const C_el = _ip_data[ip].computeElasticTangentStiffness(
            t, x_position, dt, temperature);

        auto const beta_SR =
            (1 - alpha) /
            _ip_data[ip].solid_material.getBulkModulus(t, x_position, &C_el);
        variables.grain_compressibility = beta_SR;

        variables.effective_pore_pressure = -chi_S_L * p_cap_ip;
        variables_prev.effective_pore_pressure = -chi_S_L_prev * p_cap_prev_ip;

        // Set volumetric strain rate for the general case without swelling.
        variables.volumetric_strain = Invariants::trace(_ip_data[ip].eps);
        variables_prev.volumetric_strain = Invariants::trace(B * u_prev);

        auto& phi = _ip_data[ip].porosity;
        {  // Porosity update
            variables_prev.porosity = _ip_data[ip].porosity_prev;
            phi = medium->property(MPL::PropertyType::porosity)
                      .template value<double>(variables, variables_prev,
                                              x_position, t, dt);
            variables.porosity = phi;
        }

        auto const rho_LR =
            liquid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);
        variables.density = rho_LR;
        auto const mu =
            liquid_phase.property(MPL::PropertyType::viscosity)
                .template value<double>(variables, x_position, t, dt);

        // Swelling and possibly volumetric strain rate update.
        updateSwellingStressAndVolumetricStrain<DisplacementDim>(
            _ip_data[ip], *medium, solid_phase, C_el, rho_LR, mu,
            _process_data.micro_porosity_parameters, alpha, phi, p_cap_ip,
            variables, variables_prev, x_position, t, dt);

        if (medium->hasProperty(MPL::PropertyType::transport_porosity))
        {
            if (!medium->hasProperty(MPL::PropertyType::saturation_micro))
            {
                variables_prev.transport_porosity =
                    _ip_data[ip].transport_porosity_prev;

                _ip_data[ip].transport_porosity =
                    medium->property(MPL::PropertyType::transport_porosity)
                        .template value<double>(variables, variables_prev,
                                                x_position, t, dt);
                variables.transport_porosity = _ip_data[ip].transport_porosity;
            }
        }
        else
        {
            variables.transport_porosity = phi;
        }

        // Set mechanical variables for the intrinsic permeability model
        // For stress dependent permeability.
        {
            auto const sigma_total = (_ip_data[ip].sigma_eff +
                                      alpha * chi_S_L * identity2 * p_cap_ip)
                                         .eval();
            // For stress dependent permeability.
            variables.total_stress.emplace<SymmetricTensor>(
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                    sigma_total));
        }

        variables.equivalent_plastic_strain =
            _ip_data[ip].material_state_variables->getEquivalentPlasticStrain();

        auto const K_intrinsic = MPL::formEigenTensor<DisplacementDim>(
            medium->property(MPL::PropertyType::permeability)
                .value(variables, x_position, t, dt));

        double const k_rel =
            medium->property(MPL::PropertyType::relative_permeability)
                .template value<double>(variables, x_position, t, dt);

        GlobalDimMatrixType const K_over_mu = k_rel * K_intrinsic / mu;

        auto const& sigma_eff = _ip_data[ip].sigma_eff;
        double const p_FR = -chi_S_L * p_cap_ip;
        // p_SR
        variables.solid_grain_pressure =
            p_FR - sigma_eff.dot(identity2) / (3 * (1 - phi));
        auto const rho_SR =
            solid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);
        _ip_data[ip].dry_density_solid = (1 - phi) * rho_SR;

        auto& sigma_sw = _ip_data[ip].sigma_sw;

        eps_m.noalias() =
            solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate)
                ? eps + C_el.inverse() * sigma_sw
                : eps;
        variables.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m);

        _ip_data[ip].updateConstitutiveRelation(variables, t, x_position, dt,
                                                temperature);

        auto const& b = _process_data.specific_body_force;

        // Compute the velocity
        auto const& dNdx_p = _ip_data[ip].dNdx_p;
        _ip_data[ip].v_darcy.noalias() =
            -K_over_mu * dNdx_p * p_L + rho_LR * K_over_mu * b;

        saturation_avg += S_L;
        porosity_avg += phi;
        sigma_avg += sigma_eff;
    }
    saturation_avg /= n_integration_points;
    porosity_avg /= n_integration_points;
    sigma_avg /= n_integration_points;

    (*_process_data.element_saturation)[_element.getID()] = saturation_avg;
    (*_process_data.element_porosity)[_element.getID()] = porosity_avg;

    Eigen::Map<KV>(&(*_process_data.element_stresses)[_element.getID() *
                                                      KV::RowsAtCompileTime]) =
        MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_avg);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(_element, _is_axially_symmetric, p_L,
                         *_process_data.pressure_interpolated);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
unsigned RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::getNumberOfIntegrationPoints() const
{
    return _integration_method.getNumberOfPoints();
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
typename MaterialLib::Solids::MechanicsBase<
    DisplacementDim>::MaterialStateVariables const&
RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                ShapeFunctionPressure, DisplacementDim>::
    getMaterialStateVariablesAt(unsigned integration_point) const
{
    return *_ip_data[integration_point].material_state_variables;
}
}  // namespace RichardsMechanics
}  // namespace ProcessLib
