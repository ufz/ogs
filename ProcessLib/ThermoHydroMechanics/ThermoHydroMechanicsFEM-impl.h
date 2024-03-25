/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "NumLib/NumericalStability/AdvectionMatrixAssembler.h"
#include "NumLib/NumericalStability/HydrodynamicDispersion.h"
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

    // Consistency check: if frozen liquid phase is given, then the constitutive
    // relation for ice must also be given, and vice versa.
    auto const& medium = _process_data.media_map.getMedium(_element.getID());
    if (medium->hasPhase("FrozenLiquid") !=
        (_process_data.ice_constitutive_relation != nullptr))
    {
        OGS_FATAL(
            "Frozen liquid phase is {:s} and the solid material constitutive "
            "relation for ice is {:s}. But both must be given (or both "
            "omitted).",
            medium->hasPhase("FrozenLiquid") ? "specified" : "not specified",
            _process_data.ice_constitutive_relation != nullptr
                ? "specified"
                : "not specified");
    }
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

        ip_data.N = shape_matrices_p[ip].N;
        ip_data.dNdx = shape_matrices_p[ip].dNdx;

        _secondary_data.N_u[ip] = shape_matrices_u[ip].N;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::size_t ThermoHydroMechanicsLocalAssembler<
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
            "order of the local assembler for element {:d} is different from "
            "the integration order in the initial condition.",
            _element.getID());
    }

    if (name == "sigma")
    {
        if (_process_data.initial_stress.value)
        {
            OGS_FATAL(
                "Setting initial conditions for stress from integration "
                "point data and from a parameter '{:s}' is not possible "
                "simultaneously.",
                _process_data.initial_stress.value->name);
        }

        return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
            values, _ip_data, &IpData::sigma_eff);
    }
    if (name == "epsilon_m")
    {
        return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
            values, _ip_data, &IpData::eps_m);
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

        int const element_id = _element.getID();
        DBUG(
            "The solid material of element {:d} (material ID {:d}) does not "
            "have an internal state variable called {:s}.",
            element_id, (*_process_data.material_ids)[element_id], name);
    }

    return 0;
}
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void ThermoHydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::setInitialConditionsConcrete(Eigen::VectorXd const
                                                       local_x,
                                                   double const t,
                                                   int const /*process_id*/)
{
    if (!_process_data.initial_stress.isTotalStress())
    {
        return;
    }

    // TODO: For staggered scheme, overload
    // LocalAssemblerInterface::setInitialConditions to enable local_x contains
    // the primary variables from all coupled processes.
    auto const p = local_x.template segment<pressure_size>(pressure_index);

    double const dt = 0.0;

    MPL::VariableArray vars;
    auto const& medium = _process_data.media_map.getMedium(_element.getID());

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const& N = _ip_data[ip].N;
        auto const& N_u = _ip_data[ip].N_u;
        ParameterLib::SpatialPosition const x_position{
            std::nullopt, _element.getID(), ip,
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    _element, N_u))};
        auto const alpha_b =
            medium->property(MPL::PropertyType::biot_coefficient)
                .template value<double>(vars, x_position, t, dt);

        auto& sigma_eff = _ip_data[ip].sigma_eff;
        sigma_eff.noalias() += alpha_b * N.dot(p) * Invariants::identity2;
        _ip_data[ip].sigma_eff_prev.noalias() = sigma_eff;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
ConstitutiveRelationsValues<DisplacementDim> ThermoHydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    updateConstitutiveRelations(
        Eigen::Ref<Eigen::VectorXd const> const local_x,
        Eigen::Ref<Eigen::VectorXd const> const local_x_prev,
        ParameterLib::SpatialPosition const& x_position, double const t,
        double const dt, IpData& ip_data,
        IntegrationPointDataForOutput<DisplacementDim>& ip_data_output) const
{
    assert(local_x.size() ==
           pressure_size + displacement_size + temperature_size);

    auto const [T, p, u] = localDOF(local_x);
    auto const [T_prev, p_prev, u_prev] = localDOF(local_x_prev);

    auto const& solid_material =
        MaterialLib::Solids::selectSolidConstitutiveRelation(
            _process_data.solid_materials, _process_data.material_ids,
            _element.getID());

    auto const& medium = _process_data.media_map.getMedium(_element.getID());
    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& solid_phase = medium->phase("Solid");
    auto* const frozen_liquid_phase = medium->hasPhase("FrozenLiquid")
                                          ? &medium->phase("FrozenLiquid")
                                          : nullptr;
    MaterialPropertyLib::VariableArray vars;

    auto const& N_u = ip_data.N_u;
    auto const& dNdx_u = ip_data.dNdx_u;

    auto const& N = ip_data.N;
    auto const& dNdx = ip_data.dNdx;

    auto const T_int_pt = N.dot(T);
    auto const T_prev_int_pt = N.dot(T_prev);
    double const dT_int_pt = T_int_pt - T_prev_int_pt;

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
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const eps_prev =
        B * u_prev;

    vars.temperature = T_int_pt;
    double const p_int_pt = N.dot(p);
    vars.liquid_phase_pressure = p_int_pt;

    vars.liquid_saturation = 1.0;

    auto const solid_density =
        solid_phase.property(MaterialPropertyLib::PropertyType::density)
            .template value<double>(vars, x_position, t, dt);

    auto const porosity =
        medium->property(MaterialPropertyLib::PropertyType::porosity)
            .template value<double>(vars, x_position, t, dt);
    vars.porosity = porosity;
    ip_data.porosity = porosity;

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
                vars, MaterialPropertyLib::Variable::liquid_phase_pressure,
                x_position, t, dt);

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
    crv.solid_linear_thermal_expansion_coefficient =
        MPL::formKelvinVector<DisplacementDim>(
            solid_phase
                .property(
                    MaterialPropertyLib::PropertyType::thermal_expansivity)
                .value(vars, x_position, t, dt));

    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const
        dthermal_strain =
            crv.solid_linear_thermal_expansion_coefficient * dT_int_pt;

    crv.K_pT_thermal_osmosis =
        (solid_phase.hasProperty(
             MaterialPropertyLib::PropertyType::thermal_osmosis_coefficient)
             ? MaterialPropertyLib::formEigenTensor<DisplacementDim>(
                   solid_phase
                       .property(MaterialPropertyLib::PropertyType::
                                     thermal_osmosis_coefficient)
                       .value(vars, x_position, t, dt))
             : Eigen::MatrixXd::Zero(DisplacementDim, DisplacementDim));

    GlobalDimVectorType const velocity = -crv.K_over_mu * dNdx * p -
                                         crv.K_pT_thermal_osmosis * dNdx * T +
                                         crv.K_over_mu * fluid_density * b;
    ip_data_output.velocity = velocity;

    //
    // displacement equation, displacement part
    //
    auto& eps_m = ip_data.eps_m;
    auto& eps_m_prev = ip_data.eps_m_prev;
    eps_m.noalias() = eps_m_prev + eps - eps_prev - dthermal_strain;
    vars.mechanical_strain
        .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
            eps_m);

    crv.C = ip_data.updateConstitutiveRelation(vars, t, x_position, dt,
                                               T_prev_int_pt);

    crv.rho = solid_density * (1 - porosity) + porosity * fluid_density;

    crv.beta =
        porosity * fluid_volumetric_thermal_expansion_coefficient +
        (alpha - porosity) *
            Invariants::trace(crv.solid_linear_thermal_expansion_coefficient);

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

    // Thermal conductivity is moved outside and zero matrix is passed instead
    // due to multiplication with fluid's density times specific heat capacity.
    crv.effective_thermal_conductivity.noalias() +=
        fluid_density * crv.c_f *
        NumLib::computeHydrodynamicDispersion(
            _process_data.stabilizer, _element.getID(),
            GlobalDimMatrixType::Zero(DisplacementDim, DisplacementDim),
            velocity, 0. /* phi */, 0. /* dispersivity_transversal */,
            0. /*dispersivity_longitudinal*/);

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
        MaterialPropertyLib::VariableArray vars_ice;
        double const phi_fr =
            (*medium)[MaterialPropertyLib::PropertyType::volume_fraction]
                .template value<double>(vars, x_position, t, dt);
        ip_data.phi_fr = phi_fr;

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

        double const phi_fr_prev = [&]()
        {
            MaterialPropertyLib::VariableArray vars_prev;
            vars_prev.temperature = T_prev_int_pt;
            return (*medium)[MaterialPropertyLib::PropertyType::volume_fraction]
                .template value<double>(vars_prev, x_position, t, dt);
        }();
        ip_data.phi_fr_prev = phi_fr_prev;

        // alpha_T^I
        MathLib::KelvinVector::KelvinVectorType<
            DisplacementDim> const ice_linear_thermal_expansion_coefficient =
            MPL::formKelvinVector<DisplacementDim>(
                frozen_liquid_phase
                    ->property(
                        MaterialPropertyLib::PropertyType::thermal_expansivity)
                    .value(vars, x_position, t, dt));

        MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const
            dthermal_strain_ice =
                ice_linear_thermal_expansion_coefficient * dT_int_pt;

        // alpha_{phi_I} -- linear expansion coeff. due to water-to-ice
        // transition (phase change), and related phase_change_strain term
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const
            phase_change_expansion_coefficient =
                MPL::formKelvinVector<DisplacementDim>(
                    frozen_liquid_phase
                        ->property(MaterialPropertyLib::PropertyType::
                                       phase_change_expansivity)
                        .value(vars, x_position, t, dt));

        MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const
            dphase_change_strain = phase_change_expansion_coefficient *
                                   (phi_fr - phi_fr_prev) / porosity;

        // eps0 ia a 'history variable' -- a solid matrix strain accrued
        // prior to the onset of ice forming
        auto& eps0 = ip_data.eps0;
        auto const& eps0_prev = ip_data.eps0_prev;

        // definition of eps_m_ice
        auto& eps_m_ice = ip_data.eps_m_ice;
        auto const& eps_m_ice_prev = ip_data.eps_m_ice_prev;

        eps_m_ice.noalias() = eps_m_ice_prev + eps - eps_prev -
                              (eps0 - eps0_prev) - dthermal_strain_ice -
                              dphase_change_strain;

        vars_ice.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m_ice);
        auto const C_IR = ip_data.updateConstitutiveRelationIce(
            *_process_data.ice_constitutive_relation, vars_ice, t, x_position,
            dt, T_prev_int_pt);
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

        crv.J_uu_fr = phi_fr * C_IR;

        auto const& sigma_eff_ice = ip_data.sigma_eff_ice;
        crv.r_u_fr = phi_fr * sigma_eff_ice;

        crv.J_uT_fr = phi_fr * C_IR * ice_linear_thermal_expansion_coefficient;

        crv.J_TT_fr = ((rho_fr * c_fr - fluid_density * crv.c_f) * dphi_fr_dT +
                       l_fr * rho_fr * d2phi_fr_dT2) *
                      dT_int_pt / dt;
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
                         std::vector<double> const& local_x_prev,
                         std::vector<double>& /*local_M_data*/,
                         std::vector<double>& /*local_K_data*/,
                         std::vector<double>& local_rhs_data,
                         std::vector<double>& local_Jac_data)
{
    assert(local_x.size() ==
           pressure_size + displacement_size + temperature_size);

    auto const x =
        Eigen::Map<Eigen::VectorXd const>(local_x.data(), local_x.size());
    auto const x_prev = Eigen::Map<Eigen::VectorXd const>(local_x_prev.data(),
                                                          local_x_prev.size());

    auto const [T, p, u] = localDOF(local_x);
    auto const [T_prev, p_prev, u_prev] = localDOF(local_x_prev);

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

    auto const& medium = _process_data.media_map.getMedium(_element.getID());
    bool const has_frozen_liquid_phase = medium->hasPhase("FrozenLiquid");

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    std::vector<GlobalDimVectorType> ip_flux_vector;
    double average_velocity_norm = 0.0;
    ip_flux_vector.reserve(n_integration_points);

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
            x, x_prev, x_position, t, dt, _ip_data[ip], _ip_data_output[ip]);

        auto const& w = _ip_data[ip].integration_weight;

        auto const& dNdx_u = _ip_data[ip].dNdx_u;

        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        auto const T_int_pt = N.dot(T);

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

        if (has_frozen_liquid_phase)
        {
            local_Jac
                .template block<displacement_size, displacement_size>(
                    displacement_index, displacement_index)
                .noalias() += B.transpose() * crv.J_uu_fr * B * w;

            local_rhs.template segment<displacement_size>(displacement_index)
                .noalias() -= B.transpose() * crv.r_u_fr * w;

            local_Jac
                .template block<displacement_size, temperature_size>(
                    displacement_index, temperature_index)
                .noalias() -= B.transpose() * crv.J_uT_fr * N * w;
        }

        local_Jac
            .template block<displacement_size, displacement_size>(
                displacement_index, displacement_index)
            .noalias() += B.transpose() * crv.C * B * w;

        local_Jac
            .template block<displacement_size, temperature_size>(
                displacement_index, temperature_index)
            .noalias() -= B.transpose() * crv.C *
                          crv.solid_linear_thermal_expansion_coefficient * N *
                          w;

        local_rhs.template segment<displacement_size>(displacement_index)
            .noalias() -= (B.transpose() * _ip_data[ip].sigma_eff -
                           N_u_op(N_u).transpose() * crv.rho * b) *
                          w;

        //
        // displacement equation, pressure part (K_up)
        //
        Kup.noalias() +=
            B.transpose() * crv.alpha_biot * Invariants::identity2 * N * w;

        //
        // pressure equation, pressure part (K_pp and M_pp).
        //
        laplace_p.noalias() += dNdx.transpose() * crv.K_over_mu * dNdx * w;

        storage_p.noalias() +=
            N.transpose() *
            (_ip_data[ip].porosity * crv.fluid_compressibility +
             (crv.alpha_biot - _ip_data[ip].porosity) * crv.beta_SR) *
            N * w;

        laplace_T.noalias() +=
            dNdx.transpose() * crv.K_pT_thermal_osmosis * dNdx * w;
        //
        //  RHS, pressure part
        //
        double const fluid_density = _ip_data_output[ip].fluid_density;
        local_rhs.template segment<pressure_size>(pressure_index).noalias() +=
            dNdx.transpose() * fluid_density * crv.K_over_mu * b * w;
        //
        // pressure equation, temperature part (M_pT)
        //
        storage_T.noalias() += N.transpose() * crv.beta * N * w;

        //
        // pressure equation, displacement part.
        //
        // Reusing Kup.transpose().

        //
        // temperature equation, temperature part.
        //
        KTT.noalias() +=
            dNdx.transpose() * crv.effective_thermal_conductivity * dNdx * w;

        ip_flux_vector.emplace_back(velocity * fluid_density * crv.c_f);
        average_velocity_norm += velocity.norm();

        if (has_frozen_liquid_phase)
        {
            local_Jac
                .template block<temperature_size, temperature_size>(
                    temperature_index, temperature_index)
                .noalias() -= N.transpose() * crv.J_TT_fr * N * w;
        }

        MTT.noalias() +=
            N.transpose() * crv.effective_volumetric_heat_capacity * N * w;

        //
        // temperature equation, pressure part
        //
        KTp.noalias() +=
            dNdx.transpose() * T_int_pt * crv.K_pT_thermal_osmosis * dNdx * w;

        // linearized darcy
        dKTT_dp.noalias() -= fluid_density * crv.c_f * N.transpose() *
                             (dNdx * T).transpose() * crv.K_over_mu * dNdx * w;

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
                dNdx.transpose() *
                (-T_int_pt * fluid_volumetric_thermal_expansion_coefficient *
                K_pT_thermal_osmosis / fluid_compressibility) *
                dNdx * w;

            local_rhs.template segment<temperature_size>(temperature_index)
                .noalias() +=
                dNdx.transpose() *
                (-T_int_pt * fluid_volumetric_thermal_expansion_coefficient /
                fluid_compressibility) *
                fluid_density * K_over_mu * b * w;
            MTu part for rhs and Jacobian:
                (-T_int_pt *
                Invariants::trace(solid_linear_thermal_expansion_coefficient) /
                solid_skeleton_compressibility) *
                N.transpose() * identity2.transpose() * B * w;
            KTp part for rhs and Jacobian:
                dNdx.transpose() *
                (T_int_pt * fluid_volumetric_thermal_expansion_coefficient *
                K_over_mu / fluid_compressibility) *
                dNdx * w;
        }
         */
    }

    NumLib::assembleAdvectionMatrix(
        _process_data.stabilizer, _ip_data, ip_flux_vector,
        average_velocity_norm / static_cast<double>(n_integration_points), KTT);

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
        laplace_p * p + laplace_T * T + storage_p * (p - p_prev) / dt -
        storage_T * (T - T_prev) / dt + Kup.transpose() * (u - u_prev) / dt;

    // displacement equation (f_u)
    local_rhs.template segment<displacement_size>(displacement_index)
        .noalias() += Kup * p;

    // temperature equation (f_T)
    local_rhs.template segment<temperature_size>(temperature_index).noalias() -=
        KTT * T + MTT * (T - T_prev) / dt;

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
    computeSecondaryVariableConcrete(double const /*t*/, double const /*dt*/,
                                     Eigen::VectorXd const& local_x,
                                     Eigen::VectorXd const& /*local_x_prev*/)
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
