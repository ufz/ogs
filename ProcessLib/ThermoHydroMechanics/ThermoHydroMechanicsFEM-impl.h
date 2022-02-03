/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <typeinfo>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/MPL/Utils/FormKelvinVector.h"
#include "MaterialLib/MPL/Utils/GetLiquidThermalExpansivity.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/EigenBlockMatrixView.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "ThermoHydroMechanicsFEM.h"

namespace ProcessLib
{
namespace ThermoHydroMechanics
{
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
ThermoHydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                   ShapeFunctionPressure, DisplacementDim>::
    ThermoHydroMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        ThermoHydroMechanicsProcessData<DisplacementDim>& process_data)
    : _process_data(process_data),
      _integration_method(integration_method),
      _element(e),
      _is_axially_symmetric(is_axially_symmetric)
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);
    _ip_data_output.resize(n_integration_points);
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

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        _ip_data.emplace_back(solid_material);
        auto& ip_data = _ip_data[ip];
        auto const& sm_u = shape_matrices_u[ip];
        ip_data.integration_weight =
            _integration_method.getWeightedPoint(ip).getWeight() *
            sm_u.integralMeasure * sm_u.detJ;

        ip_data.N_u = sm_u.N;
        ip_data.dNdx_u = sm_u.dNdx;

        ip_data.N_p = shape_matrices_p[ip].N;
        ip_data.dNdx_p = shape_matrices_p[ip].dNdx;

        _secondary_data.N_u[ip] = shape_matrices_u[ip].N;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::size_t ThermoHydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
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
    if (name == "epsilon_ip")
    {
        return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
            values, _ip_data, &IpData::eps);
    }

    return 0;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
ConstitutiveRelationsValues<DisplacementDim> ThermoHydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    updateConstitutiveRelations(
        Eigen::Ref<Eigen::VectorXd const> const local_x,
        Eigen::Ref<Eigen::VectorXd const> const local_xdot,
        ParameterLib::SpatialPosition const& x_position, double const t,
        double const dt, IpData& ip_data,
        IntegrationPointDataForOutput<DisplacementDim>& ip_data_output) const
{
    assert(local_x.size() ==
           pressure_size + displacement_size + temperature_size);

    auto const [T, p, u] = localDOF(local_x);
    auto const [T_dot, p_dot, u_dot] = localDOF(local_xdot);

    auto const& solid_material =
        MaterialLib::Solids::selectSolidConstitutiveRelation(
            _process_data.solid_materials, _process_data.material_ids,
            _element.getID());

    auto const& medium = _process_data.media_map->getMedium(_element.getID());
    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& solid_phase = medium->phase("Solid");
    auto* const frozen_liquid_phase = medium->hasPhase("FrozenLiquid")
                                          ? &medium->phase("FrozenLiquid")
                                          : nullptr;
    MaterialPropertyLib::VariableArray vars;

    auto const& N_u = ip_data.N_u;
    auto const& dNdx_u = ip_data.dNdx_u;

    auto const& N_p = ip_data.N_p;
    auto const& dNdx_p = ip_data.dNdx_p;

    // same shape function for pressure and temperature since they have the
    // same order
    auto const& N_T = N_p;
    auto const& dNdx_T = dNdx_p;
    auto const T_int_pt = N_T.dot(T);
    double const dT_int_pt = N_T.dot(T_dot) * dt;

    auto const x_coord =
        NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                       ShapeMatricesTypeDisplacement>(_element,
                                                                      N_u);
    auto const B =
        LinearBMatrix::computeBMatrix<DisplacementDim,
                                      ShapeFunctionDisplacement::NPOINTS,
                                      typename BMatricesType::BMatrixType>(
            dNdx_u, N_u, x_coord, _is_axially_symmetric);

    ConstitutiveRelationsValues<DisplacementDim> crv;

    auto& eps = ip_data.eps;
    eps.noalias() = B * u;

    vars.temperature = T_int_pt;
    double const p_int_pt = N_p.dot(p);
    vars.phase_pressure = p_int_pt;

    vars.liquid_saturation = 1.0;

    auto const solid_density =
        solid_phase.property(MaterialPropertyLib::PropertyType::density)
            .template value<double>(vars, x_position, t, dt);

    auto const porosity =
        medium->property(MaterialPropertyLib::PropertyType::porosity)
            .template value<double>(vars, x_position, t, dt);
    vars.porosity = porosity;
    crv.porosity = porosity;

    crv.alpha_biot =
        medium->property(MaterialPropertyLib::PropertyType::biot_coefficient)
            .template value<double>(vars, x_position, t, dt);
    auto const& alpha = crv.alpha_biot;

    auto const C_el = ip_data.computeElasticTangentStiffness(
        t, x_position, dt, static_cast<double>(T_int_pt));
    auto const solid_skeleton_compressibility =
        1 / solid_material.getBulkModulus(t, x_position, &C_el);

    crv.beta_SR = (1 - alpha) * solid_skeleton_compressibility;

    // Set mechanical variables for the intrinsic permeability model
    // For stress dependent permeability.
    {
        auto const& identity2 = Invariants::identity2;
        auto const sigma_total =
            (ip_data.sigma_eff - alpha * p_int_pt * identity2).eval();
        vars.total_stress.emplace<SymmetricTensor>(
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_total));
    }
    // For strain dependent permeability
    vars.volumetric_strain = Invariants::trace(ip_data.eps);
    vars.equivalent_plastic_strain =
        ip_data.material_state_variables->getEquivalentPlasticStrain();

    auto const intrinsic_permeability =
        MaterialPropertyLib::formEigenTensor<DisplacementDim>(
            medium->property(MaterialPropertyLib::PropertyType::permeability)
                .value(vars, x_position, t, dt));

    auto const fluid_density =
        liquid_phase.property(MaterialPropertyLib::PropertyType::density)
            .template value<double>(vars, x_position, t, dt);
    ip_data_output.fluid_density = fluid_density;
    vars.density = fluid_density;

    auto const drho_dp =
        liquid_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(
                vars, MaterialPropertyLib::Variable::phase_pressure, x_position,
                t, dt);

    crv.fluid_compressibility = 1 / fluid_density * drho_dp;

    double const fluid_volumetric_thermal_expansion_coefficient =
        MaterialPropertyLib::getLiquidThermalExpansivity(
            liquid_phase, vars, fluid_density, x_position, t, dt);

    // Use the viscosity model to compute the viscosity
    ip_data_output.viscosity =
        liquid_phase.property(MaterialPropertyLib::PropertyType::viscosity)
            .template value<double>(vars, x_position, t, dt);
    crv.K_over_mu = intrinsic_permeability / ip_data_output.viscosity;

    auto const& b = _process_data.specific_body_force;

    // Consider also anisotropic thermal expansion.
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const
        solid_linear_thermal_expansion_coefficient =
            MPL::formKelvinVector<DisplacementDim>(
                solid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_expansivity)
                    .value(vars, x_position, t, dt));

    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const
        dthermal_strain =
            solid_linear_thermal_expansion_coefficient * dT_int_pt;

    crv.K_pT_thermal_osmosis =
        (solid_phase.hasProperty(
             MaterialPropertyLib::PropertyType::thermal_osmosis_coefficient)
             ? MaterialPropertyLib::formEigenTensor<DisplacementDim>(
                   solid_phase
                       .property(MaterialPropertyLib::PropertyType::
                                     thermal_osmosis_coefficient)
                       .value(vars, x_position, t, dt))
             : Eigen::MatrixXd::Zero(DisplacementDim, DisplacementDim));

    GlobalDimVectorType const velocity = -crv.K_over_mu * dNdx_p * p -
                                         crv.K_pT_thermal_osmosis * dNdx_T * T +
                                         crv.K_over_mu * fluid_density * b;
    ip_data_output.velocity = velocity;

    //
    // displacement equation, displacement part
    //
    auto& eps_prev = ip_data.eps_prev;
    auto& eps_m = ip_data.eps_m;
    auto& eps_m_prev = ip_data.eps_m_prev;
    eps_m.noalias() = eps_m_prev + eps - eps_prev - dthermal_strain;
    vars.mechanical_strain
        .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
            eps_m);

    crv.C = ip_data.updateConstitutiveRelation(vars, t, x_position, dt,
                                               T_int_pt - dT_int_pt);

    crv.rho = solid_density * (1 - porosity) + porosity * fluid_density;

    crv.beta =
        porosity * fluid_volumetric_thermal_expansion_coefficient +
        (alpha - porosity) *
            Invariants::trace(solid_linear_thermal_expansion_coefficient);

    //
    // pressure equation, displacement part.
    //
    // Reusing Kup.transpose().

    //
    // temperature equation, temperature part.
    //
    crv.c_f =
        liquid_phase
            .property(MaterialPropertyLib::PropertyType::specific_heat_capacity)
            .template value<double>(vars, x_position, t, dt);
    crv.effective_thermal_conductivity =
        MaterialPropertyLib::formEigenTensor<DisplacementDim>(
            medium
                ->property(
                    MaterialPropertyLib::PropertyType::thermal_conductivity)
                .value(vars, x_position, t, dt));

    if (_process_data.stabilizer)
    {
        GlobalDimMatrixType const& I(
            GlobalDimMatrixType::Identity(DisplacementDim, DisplacementDim));
        crv.effective_thermal_conductivity.noalias() +=
            fluid_density * crv.c_f *
            _process_data.stabilizer->computeArtificialDiffusion(
                _element.getID(), velocity.norm()) *
            I;
    }

    double const c_s =
        solid_phase
            .property(MaterialPropertyLib::PropertyType::specific_heat_capacity)
            .template value<double>(vars, x_position, t, dt);

    // Also modified by freezing terms.
    crv.effective_volumetric_heat_capacity =
        porosity * fluid_density * crv.c_f +
        (1.0 - porosity) * solid_density * c_s;

    if (frozen_liquid_phase)
    {
        double const phi_fr =
            (*medium)[MaterialPropertyLib::PropertyType::volume_fraction]
                .template value<double>(vars, x_position, t, dt);

        auto const frozen_liquid_value =
            [&](MaterialPropertyLib::PropertyType const p)
        {
            return (*frozen_liquid_phase)[p].template value<double>(
                vars, x_position, t, dt);
        };

        double const rho_fr =
            frozen_liquid_value(MaterialPropertyLib::PropertyType::density);

        double const c_fr = frozen_liquid_value(
            MaterialPropertyLib::PropertyType::specific_heat_capacity);

        double const l_fr = frozen_liquid_value(
            MaterialPropertyLib::PropertyType::specific_latent_heat);

        double const dphi_fr_dT =
            (*medium)[MaterialPropertyLib::PropertyType::volume_fraction]
                .template dValue<double>(
                    vars, MaterialPropertyLib::Variable::temperature,
                    x_position, t, dt);

        crv.effective_volumetric_heat_capacity +=
            -phi_fr * fluid_density * crv.c_f + phi_fr * rho_fr * c_fr -
            l_fr * rho_fr * dphi_fr_dT;

        // part of dMTT_dT derivative for freezing
        double const d2phi_fr_dT2 =
            (*medium)[MaterialPropertyLib::PropertyType::volume_fraction]
                .template d2Value<double>(
                    vars, MaterialPropertyLib::Variable::temperature,
                    MaterialPropertyLib::Variable::temperature, x_position, t,
                    dt);

        crv.J_TT_fr = ((rho_fr * c_fr - fluid_density * crv.c_f) * dphi_fr_dT +
                       l_fr * rho_fr * d2phi_fr_dT2) *
                      N_T.dot(T_dot);
    }
    return crv;
}

// Assembles the local Jacobian matrix. So far, the linearisation of HT part is
// not considered as that in HT process.
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void ThermoHydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    assembleWithJacobian(double const t, double const dt,
                         std::vector<double> const& local_x,
                         std::vector<double> const& local_xdot,
                         std::vector<double>& /*local_M_data*/,
                         std::vector<double>& /*local_K_data*/,
                         std::vector<double>& local_rhs_data,
                         std::vector<double>& local_Jac_data)
{
    assert(local_x.size() ==
           pressure_size + displacement_size + temperature_size);

    auto const x =
        Eigen::Map<Eigen::VectorXd const>(local_x.data(), local_x.size());
    auto const xdot =
        Eigen::Map<Eigen::VectorXd const>(local_xdot.data(), local_xdot.size());

    auto const [T, p, u] = localDOF(local_x);
    auto const [T_dot, p_dot, u_dot] = localDOF(local_xdot);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            temperature_size + displacement_size + pressure_size,
            temperature_size + displacement_size + pressure_size>>(
        local_Jac_data, displacement_size + pressure_size + temperature_size,
        displacement_size + pressure_size + temperature_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size + pressure_size + temperature_size>>(
        local_rhs_data, displacement_size + pressure_size + temperature_size);

    typename ShapeMatricesTypePressure::NodalMatrixType MTT;
    MTT.setZero(temperature_size, temperature_size);

    typename ShapeMatricesTypePressure::NodalMatrixType KTT;
    KTT.setZero(temperature_size, temperature_size);

    typename ShapeMatricesTypePressure::NodalMatrixType K_TT_advection;
    K_TT_advection.setZero(temperature_size, temperature_size);

    typename ShapeMatricesTypePressure::NodalMatrixType KTp;
    KTp.setZero(temperature_size, pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType dKTT_dp;
    dKTT_dp.setZero(temperature_size, pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType laplace_p;
    laplace_p.setZero(pressure_size, pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType laplace_T;
    laplace_T.setZero(pressure_size, temperature_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage_p;
    storage_p.setZero(pressure_size, pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage_T;
    storage_T.setZero(pressure_size, temperature_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        displacement_size, pressure_size>
        Kup;
    Kup.setZero(displacement_size, pressure_size);

    auto const& medium = _process_data.media_map->getMedium(_element.getID());
    bool const has_frozen_liquid_phase = medium->hasPhase("FrozenLiquid");

    typename ShapeMatricesTypePressure::NodalVectorType node_flux_q;
    node_flux_q.setZero(temperature_size);

    bool const apply_full_upwind =
        _process_data.stabilizer &&
        (typeid(*_process_data.stabilizer) == typeid(NumLib::FullUpwind));

    double max_velocity_magnitude = 0.;

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& N_u = _ip_data[ip].N_u;
        ParameterLib::SpatialPosition const x_position{
            std::nullopt, _element.getID(), ip,
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    _element, N_u))};

        auto const crv = updateConstitutiveRelations(
            x, xdot, x_position, t, dt, _ip_data[ip], _ip_data_output[ip]);

        auto const& w = _ip_data[ip].integration_weight;

        auto const& dNdx_u = _ip_data[ip].dNdx_u;

        auto const& N_p = _ip_data[ip].N_p;
        auto const& dNdx_p = _ip_data[ip].dNdx_p;

        // same shape function for pressure and temperature since they have the
        // same order
        auto const& N_T = N_p;
        auto const& dNdx_T = dNdx_p;
        auto const T_int_pt = N_T.dot(T);

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                _element, N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        auto const& b = _process_data.specific_body_force;
        auto const velocity = _ip_data_output[ip].velocity;

        //
        // displacement equation, displacement part
        //
        local_Jac
            .template block<displacement_size, displacement_size>(
                displacement_index, displacement_index)
            .noalias() += B.transpose() * crv.C * B * w;

        local_rhs.template segment<displacement_size>(displacement_index)
            .noalias() -= (B.transpose() * _ip_data[ip].sigma_eff -
                           N_u_op(N_u).transpose() * crv.rho * b) *
                          w;

        //
        // displacement equation, pressure part (K_up)
        //
        Kup.noalias() +=
            B.transpose() * crv.alpha_biot * Invariants::identity2 * N_p * w;

        //
        // pressure equation, pressure part (K_pp and M_pp).
        //
        laplace_p.noalias() += dNdx_p.transpose() * crv.K_over_mu * dNdx_p * w;

        storage_p.noalias() += N_p.transpose() *
                               (crv.porosity * crv.fluid_compressibility +
                                (crv.alpha_biot - crv.porosity) * crv.beta_SR) *
                               N_p * w;

        laplace_T.noalias() +=
            dNdx_p.transpose() * crv.K_pT_thermal_osmosis * dNdx_T * w;
        //
        //  RHS, pressure part
        //
        double const fluid_density = _ip_data_output[ip].fluid_density;
        local_rhs.template segment<pressure_size>(pressure_index).noalias() +=
            dNdx_p.transpose() * fluid_density * crv.K_over_mu * b * w;
        //
        // pressure equation, temperature part (M_pT)
        //
        storage_T.noalias() += N_T.transpose() * crv.beta * N_T * w;

        //
        // pressure equation, displacement part.
        //
        // Reusing Kup.transpose().

        //
        // temperature equation, temperature part.
        //
        KTT.noalias() += dNdx_T.transpose() *
                         crv.effective_thermal_conductivity * dNdx_T * w;
        K_TT_advection.noalias() += N_T.transpose() * velocity.transpose() *
                                    dNdx_T * fluid_density * crv.c_f * w;

        if (apply_full_upwind)
        {
            node_flux_q.noalias() -=
                fluid_density * crv.c_f * velocity.transpose() * dNdx_T * w;
            max_velocity_magnitude =
                std::max(max_velocity_magnitude, velocity.norm());
        }

        if (has_frozen_liquid_phase)
        {
            local_Jac
                .template block<temperature_size, temperature_size>(
                    temperature_index, temperature_index)
                .noalias() -= N_T.transpose() * crv.J_TT_fr * N_T * w;
        }

        MTT.noalias() +=
            N_T.transpose() * crv.effective_volumetric_heat_capacity * N_T * w;

        //
        // temperature equation, pressure part
        //
        KTp.noalias() += dNdx_T.transpose() * T_int_pt *
                         crv.K_pT_thermal_osmosis * dNdx_p * w;

        // linearized darcy
        dKTT_dp.noalias() -= fluid_density * crv.c_f * N_T.transpose() *
                             (dNdx_T * T).transpose() * crv.K_over_mu * dNdx_p *
                             w;

        /* TODO (Joerg) Temperature changes due to thermal dilatation of the
         * fluid, which are usually discarded as being very small.
         * Zhou et al. (10.1016/S0020-7683(98)00089-4) states that:
         * "Biot (1956) neglected this term and it is included here for
         * completeness"
         * Keeping the code here in the case these are needed for the named
         * effects in the future.
        if (fluid_compressibility != 0)
        {
            auto const C_el = _ip_data[ip].computeElasticTangentStiffness(
                t, x_position, dt, static_cast<double>(T_int_pt));
            auto const solid_skeleton_compressibility =
                1 / solid_material.getBulkModulus(t, x_position, &C_el);
            double const fluid_volumetric_thermal_expansion_coefficient =
                MaterialPropertyLib::getLiquidThermalExpansivity(
                    liquid_phase, vars, fluid_density, x_position, t, dt);

            KTT.noalias() +=
                dNdx_T.transpose() *
                (-T_int_pt * fluid_volumetric_thermal_expansion_coefficient *
                K_pT_thermal_osmosis / fluid_compressibility) *
                dNdx_T * w;

            local_rhs.template segment<temperature_size>(temperature_index)
                .noalias() +=
                dNdx_T.transpose() *
                (-T_int_pt * fluid_volumetric_thermal_expansion_coefficient /
                fluid_compressibility) *
                fluid_density * K_over_mu * b * w;
            MTu part for rhs and Jacobian:
                (-T_int_pt *
                Invariants::trace(solid_linear_thermal_expansion_coefficient) /
                solid_skeleton_compressibility) *
                N_T.transpose() * identity2.transpose() * B * w;
            KTp part for rhs and Jacobian:
                dNdx_T.transpose() *
                (T_int_pt * fluid_volumetric_thermal_expansion_coefficient *
                K_over_mu / fluid_compressibility) *
                dNdx_p * w;
        }
         */
    }

    if (apply_full_upwind &&
        max_velocity_magnitude > _process_data.stabilizer->getCutoffVelocity())
    {
        NumLib::applyFullUpwind(node_flux_q, KTT);
    }
    else
    {
        KTT.noalias() += K_TT_advection;
    }

    // temperature equation, temperature part
    local_Jac
        .template block<temperature_size, temperature_size>(temperature_index,
                                                            temperature_index)
        .noalias() += KTT + MTT / dt;

    // temperature equation, pressure part
    local_Jac
        .template block<temperature_size, pressure_size>(temperature_index,
                                                         pressure_index)
        .noalias() += KTp + dKTT_dp;

    // displacement equation, pressure part
    local_Jac
        .template block<displacement_size, pressure_size>(displacement_index,
                                                          pressure_index)
        .noalias() -= Kup;

    // pressure equation, temperature part.
    local_Jac
        .template block<pressure_size, temperature_size>(pressure_index,
                                                         temperature_index)
        .noalias() -= storage_T / dt - laplace_T;

    // pressure equation, pressure part.
    local_Jac
        .template block<pressure_size, pressure_size>(pressure_index,
                                                      pressure_index)
        .noalias() += laplace_p + storage_p / dt;

    // pressure equation, displacement part.
    local_Jac
        .template block<pressure_size, displacement_size>(pressure_index,
                                                          displacement_index)
        .noalias() += Kup.transpose() / dt;

    // pressure equation (f_p)
    local_rhs.template segment<pressure_size>(pressure_index).noalias() -=
        laplace_p * p + laplace_T * T + storage_p * p_dot - storage_T * T_dot +
        Kup.transpose() * u_dot;

    // displacement equation (f_u)
    local_rhs.template segment<displacement_size>(displacement_index)
        .noalias() += Kup * p;

    // temperature equation (f_T)
    local_rhs.template segment<temperature_size>(temperature_index).noalias() -=
        KTT * T + MTT * T_dot;

    local_rhs.template segment<temperature_size>(temperature_index).noalias() -=
        KTp * p;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> const& ThermoHydroMechanicsLocalAssembler<
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
        cache_matrix.col(ip).noalias() = _ip_data_output[ip].velocity;
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          int DisplacementDim>
std::vector<double> const& ThermoHydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunction, DisplacementDim>::
    getIntPtFluidDensity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(
        _ip_data_output,
        &IntegrationPointDataForOutput<DisplacementDim>::fluid_density, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          int DisplacementDim>
std::vector<double> const& ThermoHydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunction, DisplacementDim>::
    getIntPtViscosity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(
        _ip_data_output,
        &IntegrationPointDataForOutput<DisplacementDim>::viscosity, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void ThermoHydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    computeSecondaryVariableConcrete(double const t, double const dt,
                                     Eigen::VectorXd const& local_x,
                                     Eigen::VectorXd const& local_x_dot)
{
    auto const p = local_x.template segment<pressure_size>(pressure_index);
    auto const T =
        local_x.template segment<temperature_size>(temperature_index);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    double fluid_density_avg = 0;
    double viscosity_avg = 0;

    using KV = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    KV sigma_avg = KV::Zero();

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto& ip_data = _ip_data[ip];
        auto const& N_u = ip_data.N_u;

        ParameterLib::SpatialPosition const x_position{
            std::nullopt, _element.getID(), ip,
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    _element, N_u))};

        updateConstitutiveRelations(local_x, local_x_dot, x_position, t, dt,
                                    _ip_data[ip], _ip_data_output[ip]);

        fluid_density_avg += _ip_data_output[ip].fluid_density;
        viscosity_avg += _ip_data_output[ip].viscosity;
        sigma_avg += ip_data.sigma_eff;
    }

    fluid_density_avg /= n_integration_points;
    viscosity_avg /= n_integration_points;
    sigma_avg /= n_integration_points;

    (*_process_data.element_fluid_density)[_element.getID()] =
        fluid_density_avg;
    (*_process_data.element_viscosity)[_element.getID()] = viscosity_avg;

    Eigen::Map<KV>(&(*_process_data.element_stresses)[_element.getID() *
                                                      KV::RowsAtCompileTime]) =
        MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_avg);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(_element, _is_axially_symmetric, p,
                         *_process_data.pressure_interpolated);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(_element, _is_axially_symmetric, T,
                         *_process_data.temperature_interpolated);
}
}  // namespace ThermoHydroMechanics
}  // namespace ProcessLib
