// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <limits>
#include <typeinfo>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Phase.h"
#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/MPL/Utils/FormKelvinVector.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/EigenBlockMatrixView.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/Fem/Interpolation.h"
#include "NumLib/NumericalStability/AdvectionMatrixAssembler.h"
#include "NumLib/NumericalStability/HydrodynamicDispersion.h"
#include "ProcessLib/Common/ThermoOsmosis/ThermoOsmoticCoefficient.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "ThermoHydroMechanicsFEM.h"

namespace ProcessLib
{
namespace ThermoHydroMechanics
{

namespace
{
/// Pressure storage coefficient of the mass-balance accumulation term,
/// \f[ \phi\,c_f + (\alpha - \phi)\,\beta_{SR}
///    + \Delta S_p^{\mathrm{fr}}, \f]
/// \param porosity \f$\phi\f$, the porosity.
/// \param fluid_compressibility \f$c_f\f$, the fluid compressibility.
/// \param alpha_biot \f$\alpha\f$, the Biot coefficient.
/// \param beta_SR \f$\beta_{SR}\f$.
/// \param frozen_liquid_contribution \f$\Delta S_p^{\mathrm{fr}}\f$, the
/// storage contribution of a frozen liquid phase
/// (HydraulicPropertyValues::storage_p_fr); zero without a frozen liquid
/// phase.
double massBalancePressureStorageCoefficient(
    double const porosity, double const fluid_compressibility,
    double const alpha_biot, double const beta_SR,
    double const frozen_liquid_contribution)
{
    return porosity * fluid_compressibility +
           (alpha_biot - porosity) * beta_SR + frozen_liquid_contribution;
}

/// Temperature storage coefficient of the mass-balance accumulation term,
/// \f$\beta + \Delta S_T^{\mathrm{fr}}\f$.
/// \param beta \f$\beta\f$, the thermal storage factor
/// (HydraulicPropertyValues::beta).
/// \param frozen_liquid_contribution \f$\Delta S_T^{\mathrm{fr}}\f$, the
/// storage contribution of a frozen liquid phase
/// (HydraulicPropertyValues::storage_T_fr); zero without a frozen liquid
/// phase.
double massBalanceTemperatureStorageCoefficient(
    double const beta, double const frozen_liquid_contribution)
{
    return beta + frozen_liquid_contribution;
}

/// Rate bracket of the mass-balance accumulation term,
/// \f$S_p\,\dot p - S_T\,\dot T + \dot e\f$. Charged by the mass-balance
/// residual (via storage_p, storage_T and the volumetric strain rate term)
/// and reused for the fluid-density derivative of the residual, so the two
/// cannot diverge.
/// \param storage_p_coefficient \f$S_p\f$, the pressure storage coefficient.
/// \param p_dot \f$\dot p\f$, the pressure rate.
/// \param temperature_storage_coefficient \f$S_T\f$, the temperature storage
/// coefficient.
/// \param T_dot \f$\dot T\f$, the temperature rate.
/// \param volumetric_strain_rate_term \f$\dot e\f$, the volumetric strain
/// rate term (scheme specific).
double massBalanceStorageRate(double const storage_p_coefficient,
                              double const p_dot,
                              double const temperature_storage_coefficient,
                              double const T_dot,
                              double const volumetric_strain_rate_term)
{
    return storage_p_coefficient * p_dot -
           temperature_storage_coefficient * T_dot +
           volumetric_strain_rate_term;
}

/// Returns the derivative of the mass-balance residual with respect to the
/// fluid density, \f$\partial r_p/\partial\rho_{\mathrm L}\f$: the
/// derivative of the storage rate bracket and the derivative of the viscous
/// Darcy flux
/// \f$\rho_{\mathrm L}\,K/\mu\,k_{\mathrm{rel}}
/// (\nabla p + \rho_{\mathrm L}\boldsymbol g)\f$ with respect to the fluid
/// density. The caller multiplies it by
/// \f$\partial\rho_{\mathrm L}/\partial x\f$ for each differentiation
/// direction \f$x\f$ (pressure or temperature).
///
/// (TODO (WW): the contributions of the second fluid-density derivatives,
/// \f$\partial^2\rho_{\mathrm L}/\partial x\,\partial x'\f$, are currently
/// not assembled.)
/// \param ip_data the integration point data (shape functions and weight).
/// \param storage_rate the storage rate bracket (see
/// massBalanceStorageRate()).
/// \param velocity the Darcy velocity.
/// \param buoyancy_flux the buoyancy flux term.
template <typename IpData, typename Velocity, typename BuoyancyFlux>
auto massBalanceJacobianWrtFluidDensity(IpData const& ip_data,
                                        double const storage_rate,
                                        Velocity const& velocity,
                                        BuoyancyFlux const& buoyancy_flux)
{
    auto const& N = ip_data.N;
    auto const& dNdx = ip_data.dNdx;
    double const w = ip_data.integration_weight;

    // Derivative of the viscous Darcy flux with respect to the fluid density.
    auto const dviscous_flux_drho = (-velocity - buoyancy_flux);

    return (N.transpose() * N * (storage_rate * w) +
            dNdx.transpose() * dviscous_flux_drho * N * w)
        .eval();
}
}  // namespace

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
    if (medium->hasPhase(MaterialPropertyLib::PhaseName::FrozenLiquid) !=
        (_process_data.ice_constitutive_relation != nullptr))
    {
        OGS_FATAL(
            "Frozen liquid phase is {:s} and the solid material constitutive "
            "relation for ice is {:s}. But both must be given (or both "
            "omitted).",
            medium->hasPhase(MaterialPropertyLib::PhaseName::FrozenLiquid)
                ? "specified"
                : "not specified",
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
auto ThermoHydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::getMediumPhases() const -> MediumPhases
{
    auto const* const medium =
        _process_data.media_map.getMedium(_element.getID());
    return {medium,
            medium->phase(MaterialPropertyLib::PhaseName::AqueousLiquid),
            medium->phase(MaterialPropertyLib::PhaseName::Solid),
            medium->hasPhase(MaterialPropertyLib::PhaseName::FrozenLiquid)
                ? &medium->phase(MaterialPropertyLib::PhaseName::FrozenLiquid)
                : nullptr};
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
    // TODO: For staggered scheme, overload
    // LocalAssemblerInterface::setInitialConditions to enable local_x contains
    // the primary variables from all coupled processes.
    auto const [T, p, u] = localDOF(local_x);

    constexpr double dt = std::numeric_limits<double>::quiet_NaN();

    MPL::VariableArray vars;
    [[maybe_unused]] auto const& [medium, liquid_phase, solid_phase,
                                  frozen_liquid_phase] = getMediumPhases();

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto& ip_data = _ip_data[ip];
        auto const& N = ip_data.N;
        auto const& N_u = ip_data.N_u;
        ParameterLib::SpatialPosition const x_position{
            std::nullopt, _element.getID(),
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    _element, N_u))};

        auto& sigma_eff = ip_data.sigma_eff;
        if (_process_data.initial_stress.isTotalStress())
        {
            auto const alpha_b =
                medium->property(MPL::PropertyType::biot_coefficient)
                    .template value<double>(vars, x_position, t, dt);

            sigma_eff.noalias() += alpha_b * N.dot(p) * Invariants::identity2;
        }
        ip_data.sigma_eff_prev.noalias() = sigma_eff;

        vars.temperature = N.dot(T);
        if (frozen_liquid_phase)
        {
            auto const porosity =
                (*medium)[MaterialPropertyLib::PropertyType::porosity]
                    .template value<double>(vars, x_position, t, dt);
            auto const S_fr =
                (*medium)[MaterialPropertyLib::PropertyType::
                              frozen_liquid_saturation]
                    .template value<double>(vars, x_position, t, dt);
            ip_data.phi_fr = S_fr * porosity;
        }
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
CommonPropertyValues<DisplacementDim> ThermoHydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    updateCommonProperties(Eigen::Ref<Eigen::VectorXd const> const local_x,
                           Eigen::Ref<Eigen::VectorXd const> const local_x_prev,
                           SpaceTimeData const& space_time_data,
                           IpState const& ip_state) const
{
    auto const& [x_position, t, dt] = space_time_data;
    auto& ip_data = ip_state.ip_data;
    auto& ip_data_output = ip_state.ip_data_output;

    auto const [T, p, u] = localDOF(local_x);
    auto const [T_prev, p_prev, u_prev] = localDOF(local_x_prev);

    [[maybe_unused]] auto const& [medium, liquid_phase, solid_phase,
                                  frozen_liquid_phase] = getMediumPhases();

    MaterialPropertyLib::VariableArray vars;

    auto const& N_u = ip_data.N_u;
    auto const& dNdx_u = ip_data.dNdx_u;
    auto const& N = ip_data.N;

    CommonPropertyValues<DisplacementDim> cpv;

    auto const T_int_pt = N.dot(T);
    cpv.T_prev_int_pt = N.dot(T_prev);
    cpv.dT_int_pt = T_int_pt - cpv.T_prev_int_pt;
    cpv.T_dot = cpv.dT_int_pt / dt;
    vars.temperature = T_int_pt;

    double const p_int_pt = N.dot(p);
    cpv.p_prev_int_pt = N.dot(p_prev);
    cpv.p_dot = (p_int_pt - cpv.p_prev_int_pt) / dt;
    vars.liquid_phase_pressure = p_int_pt;

    // x_coord is used as the radial coordinate r for axisymmetric problems.
    // If axisymmetry is not enabled, this value is ignored.
    auto const x_coord = x_position.getCoordinates().value()[0];
    auto const B =
        LinearBMatrix::computeBMatrix<DisplacementDim,
                                      ShapeFunctionDisplacement::NPOINTS,
                                      typename BMatricesType::BMatrixType>(
            dNdx_u, N_u, x_coord, _is_axially_symmetric);

    auto& eps = ip_data.eps;
    eps.noalias() = B * u;
    cpv.eps_prev = B * u_prev;
    cpv.eps_v_dot =
        (Invariants::trace(eps) - Invariants::trace(cpv.eps_prev)) / dt;

    vars.liquid_saturation = 1.0;

    // Note: if a volumetric strain dependent porosity model is available in
    // future, set vars.volumetric_strain = Invariants::trace(eps) here, so that
    // the porosity property can use it to compute the porosity.
    auto const porosity =
        medium->property(MaterialPropertyLib::PropertyType::porosity)
            .template value<double>(vars, x_position, t, dt);
    vars.porosity = porosity;

    cpv.alpha_biot =
        medium->property(MaterialPropertyLib::PropertyType::biot_coefficient)
            .template value<double>(vars, x_position, t, dt);
    cpv.modified_biot = cpv.alpha_biot;

    // Set mechanical variables for the intrinsic permeability model
    // For stress dependent permeability.
    {
        auto const& identity2 = Invariants::identity2;
        auto const sigma_total =
            (ip_data.sigma_eff - cpv.alpha_biot * p_int_pt * identity2).eval();
        vars.total_stress.emplace<SymmetricTensor>(
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_total));
    }
    // For strain dependent permeability
    vars.volumetric_strain = Invariants::trace(eps);
    vars.equivalent_plastic_strain =
        ip_data.material_state_variables->getEquivalentPlasticStrain();

    auto const fluid_density =
        liquid_phase.property(MaterialPropertyLib::PropertyType::density)
            .template value<double>(vars, x_position, t, dt);
    ip_data_output.fluid_density = fluid_density;
    vars.density = fluid_density;

    cpv.drho_LR_dT =
        liquid_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(vars,
                                     MaterialPropertyLib::Variable::temperature,
                                     x_position, t, dt);

    if (frozen_liquid_phase)
    {
        auto const& frozen_liquid_saturation = (*medium)
            [MaterialPropertyLib::PropertyType::frozen_liquid_saturation];

        ip_data.phi_fr = frozen_liquid_saturation.template value<double>(
                             vars, x_position, t, dt) *
                         porosity;
        // Set frozen_liquid_saturation variable for relative permeability
        // calculation; frozen_liquid_saturation = fraction of pore space
        // occupied by ice
        vars.frozen_liquid_saturation = ip_data.phi_fr / porosity;

        double const rho_fr =
            (*frozen_liquid_phase)[MaterialPropertyLib::PropertyType::density]
                .template value<double>(vars, x_position, t, dt);
        ip_data_output.rho_fr = rho_fr;

        cpv.modified_biot += cpv.alpha_biot * ip_data.phi_fr / porosity *
                             (rho_fr / fluid_density - 1);

        ip_data.phi_fr_prev = [&]()
        {
            MaterialPropertyLib::VariableArray vars_prev;
            vars_prev.temperature = cpv.T_prev_int_pt;
            return frozen_liquid_saturation.template value<double>(
                       vars_prev, x_position, t, dt) *
                   porosity;
        }();

        cpv.dphi_fr_dT = frozen_liquid_saturation.template dValue<double>(
                             vars, MaterialPropertyLib::Variable::temperature,
                             x_position, t, dt) *
                         porosity;
        cpv.d2phi_fr_dT2 =
            frozen_liquid_saturation.template d2Value<double>(
                vars, MaterialPropertyLib::Variable::temperature,
                MaterialPropertyLib::Variable::temperature, x_position, t, dt) *
            porosity;

        // alpha_T^I
        cpv.ice_linear_thermal_expansion_coefficient =
            MPL::formKelvinVector<DisplacementDim>(
                frozen_liquid_phase
                    ->property(
                        MaterialPropertyLib::PropertyType::thermal_expansivity)
                    .value(vars, x_position, t, dt));
    }

    // Consider also anisotropic thermal expansion.
    cpv.solid_linear_thermal_expansion_coefficient =
        MPL::formKelvinVector<DisplacementDim>(
            solid_phase
                .property(
                    MaterialPropertyLib::PropertyType::thermal_expansivity)
                .value(vars, x_position, t, dt));
    cpv.solid_volumetric_thermal_expansion_coefficient =
        Invariants::trace(cpv.solid_linear_thermal_expansion_coefficient);

    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const
        dthermal_strain =
            cpv.solid_linear_thermal_expansion_coefficient * cpv.dT_int_pt;

    //
    // The mechanical strain is set here, before the variable array is passed
    // on in cpv.vars, so that permeability models that depend on it can read
    // it in the property evaluations later.
    //
    auto& eps_m = ip_data.eps_m;
    auto const& eps_m_prev = ip_data.eps_m_prev;
    eps_m.noalias() = eps_m_prev + eps - cpv.eps_prev - dthermal_strain;
    vars.mechanical_strain
        .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
            eps_m);

    cpv.solid_density =
        solid_phase.property(MaterialPropertyLib::PropertyType::density)
            .template value<double>(vars, x_position, t, dt);

    cpv.vars = vars;
    return cpv;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
FlowPropertyValues<DisplacementDim> ThermoHydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    updateFlowProperties(CommonPropertyValues<DisplacementDim> const& cpv,
                         Eigen::Ref<Eigen::VectorXd const> const local_x,
                         SpaceTimeData const& space_time_data,
                         IpState const& ip_state) const
{
    auto const& [x_position, t, dt] = space_time_data;
    auto& ip_data = ip_state.ip_data;
    auto& ip_data_output = ip_state.ip_data_output;

    auto const [T, p, u] = localDOF(local_x);

    [[maybe_unused]] auto const& [medium, liquid_phase, solid_phase,
                                  frozen_liquid_phase] = getMediumPhases();
    bool const has_frozen_liquid_phase = frozen_liquid_phase != nullptr;

    auto const& vars = cpv.vars;
    auto const& dNdx = ip_data.dNdx;

    FlowPropertyValues<DisplacementDim> fpv;

    auto const intrinsic_permeability =
        MaterialPropertyLib::formEigenTensor<DisplacementDim>(
            medium->property(MaterialPropertyLib::PropertyType::permeability)
                .value(vars, x_position, t, dt));

    ip_data_output.viscosity =
        liquid_phase.property(MaterialPropertyLib::PropertyType::viscosity)
            .template value<double>(vars, x_position, t, dt);

    fpv.K_over_mu = intrinsic_permeability / ip_data_output.viscosity;

    fpv.k_rel =
        has_frozen_liquid_phase
            ? liquid_phase
                  .property(
                      MaterialPropertyLib::PropertyType::relative_permeability)
                  .template value<double>(vars, x_position, t, dt)
            : 1.0;

    fpv.K_pT_thermal_osmosis =
        ProcessLib::getThermoOsmoticCoefficient<DisplacementDim>(
            *medium, vars, x_position, t, dt, intrinsic_permeability,
            ip_data_output.viscosity);

    fpv.grad_p = dNdx * p;

    auto const& b = _process_data.specific_body_force;
    fpv.buoyancy_flux =
        ip_data_output.fluid_density * fpv.k_rel * fpv.K_over_mu * b;
    ip_data_output.velocity = -fpv.k_rel * fpv.K_over_mu * fpv.grad_p -
                              fpv.K_pT_thermal_osmosis * dNdx * T +
                              fpv.buoyancy_flux;

    return fpv;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
HeatCapacityValues ThermoHydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    computeHeatCapacities(CommonPropertyValues<DisplacementDim> const& cpv,
                          SpaceTimeData const& space_time_data,
                          IpState const& ip_state) const
{
    auto const& [x_position, t, dt] = space_time_data;
    auto& ip_data = ip_state.ip_data;
    auto& ip_data_output = ip_state.ip_data_output;

    [[maybe_unused]] auto const& [medium, liquid_phase, solid_phase,
                                  frozen_liquid_phase] = getMediumPhases();

    auto const& vars = cpv.vars;
    double const porosity = vars.porosity;
    double const fluid_density = vars.density;

    HeatCapacityValues hcv;
    hcv.c_f =
        liquid_phase
            .property(MaterialPropertyLib::PropertyType::specific_heat_capacity)
            .template value<double>(vars, x_position, t, dt);
    hcv.c_s =
        solid_phase
            .property(MaterialPropertyLib::PropertyType::specific_heat_capacity)
            .template value<double>(vars, x_position, t, dt);

    // Sensible part (no latent term here).
    hcv.sensible_volumetric_heat_capacity =
        porosity * fluid_density * hcv.c_f +
        (1.0 - porosity) * cpv.solid_density * hcv.c_s;

    if (frozen_liquid_phase)
    {
        double const phi_fr = ip_data.phi_fr;
        hcv.c_fr =
            (*frozen_liquid_phase)
                [MaterialPropertyLib::PropertyType::specific_heat_capacity]
                    .template value<double>(vars, x_position, t, dt);

        // Freezing modifies only the sensible part; latent kept separate.
        hcv.sensible_volumetric_heat_capacity +=
            -phi_fr * fluid_density * hcv.c_f +
            phi_fr * ip_data_output.rho_fr * hcv.c_fr;
    }
    return hcv;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
ThermalPropertyValues<DisplacementDim> ThermoHydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    updateThermalProperties(CommonPropertyValues<DisplacementDim> const& cpv,
                            FlowPropertyValues<DisplacementDim> const& fpv,
                            SpaceTimeData const& space_time_data,
                            IpState const& ip_state) const
{
    auto const& [x_position, t, dt] = space_time_data;
    auto& ip_data = ip_state.ip_data;
    auto& ip_data_output = ip_state.ip_data_output;

    [[maybe_unused]] auto const& [medium, liquid_phase, solid_phase,
                                  frozen_liquid_phase] = getMediumPhases();

    auto const& vars = cpv.vars;
    ThermalPropertyValues<DisplacementDim> crv_T;

    double const porosity = vars.porosity;
    double const fluid_density = vars.density;

    // Compute relative permeability derivative with respect to temperature
    double dk_rel_dT = 0;
    if (frozen_liquid_phase)
    {
        // dk_rel/dT = (dk_rel/dphi_ice) * (dphi_ice/dT)
        // where dphi_ice/dT = dphi_fr_dT / porosity
        auto const dk_rel_dphi_ice =
            liquid_phase
                .property(
                    MaterialPropertyLib::PropertyType::relative_permeability)
                .template dValue<double>(
                    vars,
                    MaterialPropertyLib::Variable::frozen_liquid_saturation,
                    x_position, t, dt);
        dk_rel_dT = dk_rel_dphi_ice * cpv.dphi_fr_dT / porosity;
    }
    crv_T.dk_rel_dT = dk_rel_dT;

    crv_T.dvelocity_dT =
        -dk_rel_dT * fpv.K_over_mu * fpv.grad_p +
        // TODO(naumov): - fpv.K_pT_thermal_osmosis * dNdx * dT_dT +
        (dk_rel_dT * fluid_density + fpv.k_rel * cpv.drho_LR_dT) *
            fpv.K_over_mu * _process_data.specific_body_force;

    //
    // temperature equation, temperature part.
    //
    auto const hcv = computeHeatCapacities(cpv, space_time_data, ip_state);
    crv_T.c_f = hcv.c_f;
    crv_T.sensible_volumetric_heat_capacity =
        hcv.sensible_volumetric_heat_capacity;

    crv_T.effective_thermal_conductivity =
        MaterialPropertyLib::formEigenTensor<DisplacementDim>(
            medium
                ->property(
                    MaterialPropertyLib::PropertyType::thermal_conductivity)
                .value(vars, x_position, t, dt));

    crv_T.dlambda_eff_dT =
        MaterialPropertyLib::formEigenTensor<DisplacementDim>(
            medium
                ->property(
                    MaterialPropertyLib::PropertyType::thermal_conductivity)
                .dValue(vars, MaterialPropertyLib::Variable::temperature,
                        x_position, t, dt));

    // Thermal conductivity is moved outside and zero matrix is passed instead
    // due to multiplication with fluid's density times specific heat capacity.
    crv_T.effective_thermal_conductivity.noalias() +=
        fluid_density * crv_T.c_f *
        NumLib::computeHydrodynamicDispersion(
            _process_data.stabilizer, _element.getID(),
            GlobalDimMatrixType::Zero(DisplacementDim, DisplacementDim),
            ip_data_output.velocity, 0. /* phi */,
            0. /* dispersivity_transversal */,
            0. /*dispersivity_longitudinal*/);

    auto const drho_SR_dT =
        solid_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(vars,
                                     MaterialPropertyLib::Variable::temperature,
                                     x_position, t, dt);

    double dC_eff_dT = porosity * cpv.drho_LR_dT * crv_T.c_f +
                       (1.0 - porosity) * drho_SR_dT * hcv.c_s;

    if (frozen_liquid_phase)
    {
        double const phi_fr = ip_data.phi_fr;

        double const c_fr = hcv.c_fr;
        double const l_fr =
            (*frozen_liquid_phase)
                [MaterialPropertyLib::PropertyType::specific_latent_heat]
                    .template value<double>(vars, x_position, t, dt);

        double const rho_fr = ip_data_output.rho_fr;

        // Latent contribution L (effective = sensible - L).
        crv_T.latent_volumetric_heat_capacity = l_fr * rho_fr * cpv.dphi_fr_dT;

        // part of dMTT_dT derivative for freezing
        dC_eff_dT += -cpv.dphi_fr_dT * fluid_density * crv_T.c_f -
                     phi_fr * cpv.drho_LR_dT * crv_T.c_f +
                     cpv.dphi_fr_dT * rho_fr * c_fr -
                     l_fr * rho_fr * cpv.d2phi_fr_dT2;
    }
    crv_T.J_TT = dC_eff_dT * cpv.dT_int_pt / dt;
    return crv_T;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
double ThermoHydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    solidSkeletonCompressibility(SpaceTimeData const& space_time_data,
                                 IpData& ip_data,
                                 double const temperature) const
{
    auto const& [x_position, t, dt] = space_time_data;
    auto const& solid_material =
        MaterialLib::Solids::selectSolidConstitutiveRelation(
            _process_data.solid_materials, _process_data.material_ids,
            _element.getID());

    auto const C_el =
        ip_data.computeElasticTangentStiffness(space_time_data, temperature);
    return 1 / solid_material.getBulkModulus(t, x_position, &C_el);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
HydraulicPropertyValues ThermoHydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    updateHydraulicProperties(CommonPropertyValues<DisplacementDim> const& cpv,
                              SpaceTimeData const& space_time_data,
                              IpState const& ip_state) const
{
    auto const& [x_position, t, dt] = space_time_data;
    auto& ip_data = ip_state.ip_data;
    auto& ip_data_output = ip_state.ip_data_output;

    [[maybe_unused]] auto const& [medium, liquid_phase, solid_phase,
                                  frozen_liquid_phase] = getMediumPhases();
    bool const has_frozen_liquid_phase = frozen_liquid_phase != nullptr;

    auto const& vars = cpv.vars;

    HydraulicPropertyValues crv_p;

    crv_p.solid_skeleton_compressibility = solidSkeletonCompressibility(
        space_time_data, ip_data, static_cast<double>(vars.temperature));

    double const porosity = vars.porosity;
    double const fluid_density = vars.density;

    auto const drho_dp =
        liquid_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(
                vars, MaterialPropertyLib::Variable::liquid_phase_pressure,
                x_position, t, dt);
    crv_p.drho_L_dp = drho_dp;

    crv_p.fluid_compressibility = 1 / fluid_density * drho_dp;

    crv_p.beta_SR = (1 - cpv.alpha_biot) * crv_p.solid_skeleton_compressibility;

    double const alpha_minus_porosity = cpv.alpha_biot - porosity;

    double const fluid_volumetric_thermal_expansion_coefficient =
        (fluid_density == 0.0) ? 0.0 : -cpv.drho_LR_dT / fluid_density;
    crv_p.beta = porosity * fluid_volumetric_thermal_expansion_coefficient +
                 alpha_minus_porosity *
                     cpv.solid_volumetric_thermal_expansion_coefficient;

    if (has_frozen_liquid_phase)
    {
        double const phi_fr = ip_data.phi_fr;
        double const rho_fr = ip_data_output.rho_fr;
        double const dphi_fr_dT = cpv.dphi_fr_dT;

        double const mass_exchange =
            -dphi_fr_dT * porosity * (1. - rho_fr / fluid_density);

        double const beta_T_SI =
            porosity * Invariants::trace(
                           cpv.ice_linear_thermal_expansion_coefficient) +
            alpha_minus_porosity *
                cpv.solid_volumetric_thermal_expansion_coefficient;

        auto const C_el_ice = ip_data.computeElasticTangentStiffnessIce(
            *_process_data.ice_constitutive_relation, space_time_data,
            static_cast<double>(vars.temperature));
        double const beta_IR =
            1. / _process_data.ice_constitutive_relation->getBulkModulus(
                     t, x_position, &C_el_ice);

        double const storage_p_fr_coeff =
            (porosity * beta_IR + alpha_minus_porosity * crv_p.beta_SR) *
                rho_fr / fluid_density -
            (porosity * crv_p.fluid_compressibility +
             alpha_minus_porosity * crv_p.beta_SR);
        crv_p.storage_p_fr = phi_fr / porosity * storage_p_fr_coeff;

        crv_p.storage_T_fr =
            phi_fr / porosity *
                (beta_T_SI * rho_fr / fluid_density - crv_p.beta) -
            mass_exchange;

        double const dstorage_p_fr_coeff_dT =
            (porosity * beta_IR + alpha_minus_porosity * crv_p.beta_SR) *
            cpv.drho_LR_dT / (fluid_density * fluid_density);

        double const dmass_exchange_dT =
            -cpv.d2phi_fr_dT2 * porosity * (1. - rho_fr / fluid_density) +
            dphi_fr_dT * porosity * rho_fr * cpv.drho_LR_dT /
                (fluid_density * fluid_density);

        double const dstorage_T_fr_dT =
            dphi_fr_dT / porosity *
                (beta_T_SI * rho_fr / fluid_density - crv_p.beta) +
            phi_fr / porosity * beta_T_SI * rho_fr * cpv.drho_LR_dT /
                (fluid_density * fluid_density) -
            dmass_exchange_dT;

        crv_p.J_pT_fr = (dphi_fr_dT * storage_p_fr_coeff +
                         phi_fr * dstorage_p_fr_coeff_dT) /
                        porosity * cpv.p_dot;
        crv_p.J_pT_fr += dstorage_T_fr_dT * cpv.T_dot;
    }
    return crv_p;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
MechanicalPropertyValues<DisplacementDim> ThermoHydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    updateMechanicalProperties(CommonPropertyValues<DisplacementDim> const& cpv,
                               SpaceTimeData const& space_time_data,
                               IpState const& ip_state) const
{
    auto const& [x_position, t, dt] = space_time_data;
    auto& ip_data = ip_state.ip_data;
    auto& ip_data_output = ip_state.ip_data_output;

    [[maybe_unused]] auto const& [medium, liquid_phase, solid_phase,
                                  frozen_liquid_phase] = getMediumPhases();
    auto const& vars = cpv.vars;

    MechanicalPropertyValues<DisplacementDim> crv_u;

    double const porosity = vars.porosity;
    double const fluid_density = vars.density;

    crv_u.C = ip_data.updateConstitutiveRelation(vars, space_time_data,
                                                 cpv.T_prev_int_pt);

    crv_u.rho = cpv.solid_density * (1 - porosity) + porosity * fluid_density;

    if (frozen_liquid_phase)
    {
        double const phi_fr = ip_data.phi_fr;
        double const rho_fr = ip_data_output.rho_fr;

        crv_u.rho += phi_fr * rho_fr - phi_fr * fluid_density;

        MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const
            dthermal_strain_ice =
                cpv.ice_linear_thermal_expansion_coefficient * cpv.dT_int_pt;

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
                                   (phi_fr - ip_data.phi_fr_prev) / porosity;

        // eps0 ia a 'history variable' -- a solid matrix strain accrued
        // prior to the onset of ice forming
        auto const& eps0 = ip_data.eps0;
        auto const& eps0_prev = ip_data.eps0_prev;

        // definition of eps_m_ice
        auto& eps_m_ice = ip_data.eps_m_ice;
        auto const& eps_m_ice_prev = ip_data.eps_m_ice_prev;

        eps_m_ice.noalias() = eps_m_ice_prev + ip_data.eps - cpv.eps_prev -
                              (eps0 - eps0_prev) - dthermal_strain_ice -
                              dphase_change_strain;

        MaterialPropertyLib::VariableArray vars_ice;
        vars_ice.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m_ice);
        auto const C_IR = ip_data.updateConstitutiveRelationIce(
            *_process_data.ice_constitutive_relation, vars_ice, space_time_data,
            cpv.T_prev_int_pt);

        crv_u.J_uu_fr = phi_fr * C_IR;
        crv_u.J_uT_fr =
            phi_fr * C_IR * cpv.ice_linear_thermal_expansion_coefficient;
        crv_u.r_u_fr = phi_fr * ip_data.sigma_eff_ice;
    }

    return crv_u;
}

// Assembles the local Jacobian matrix. So far, the linearisation of HT part
// is not considered as that in HT process.
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void ThermoHydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::assembleWithJacobian(double const t, double const dt,
                                           std::vector<double> const& local_x,
                                           std::vector<double> const&
                                               local_x_prev,
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

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    EnergyBalanceTerms energy_balance(n_integration_points);

    typename ShapeMatricesTypePressure::NodalMatrixType dKTT_dp_T;
    dKTT_dp_T.setZero(temperature_size, pressure_size);

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

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        pressure_size, displacement_size>
        Kpu;
    if (!_process_data.is_volume_balance_equation_type)
    {
        Kpu.setZero(pressure_size, displacement_size);
    }

    [[maybe_unused]] auto const& [medium, liquid_phase, solid_phase,
                                  frozen_liquid_phase] = getMediumPhases();
    bool const has_frozen_liquid_phase = frozen_liquid_phase != nullptr;

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto& ip_data = _ip_data[ip];
        auto const& N_u = ip_data.N_u;
        auto const space_time_data = makeSpaceTimeData(ip_data.N_u, t, dt);
        IpState const ip_state{ip_data, _ip_data_output[ip]};

        // Get common properties
        auto const cpv =
            updateCommonProperties(x, x_prev, space_time_data, ip_state);
        auto const fpv =
            updateFlowProperties(cpv, x, space_time_data, ip_state);
        auto const crv_T =
            updateThermalProperties(cpv, fpv, space_time_data, ip_state);
        auto const crv_p =
            updateHydraulicProperties(cpv, space_time_data, ip_state);
        auto const crv_u =
            updateMechanicalProperties(cpv, space_time_data, ip_state);

        auto const& w = ip_data.integration_weight;

        auto const& dNdx_u = ip_data.dNdx_u;

        auto const& N = ip_data.N;
        auto const& dNdx = ip_data.dNdx;

        auto const x_coord = space_time_data.x.getCoordinates().value()[0];
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        auto const& b = _process_data.specific_body_force;
        auto const velocity = _ip_data_output[ip].velocity;

        //
        // displacement equation, displacement part (K_uu and f_u).
        //
        addMomentumBalanceStiffnessAndResidual(
            local_Jac.template block<displacement_size, displacement_size>(
                displacement_index, displacement_index),
            local_rhs.template segment<displacement_size>(displacement_index),
            ip_data, crv_u, B, ip_data.sigma_eff, has_frozen_liquid_phase);

        auto const uT_coeff =
            has_frozen_liquid_phase
                ? (crv_u.J_uT_fr +
                   crv_u.C * cpv.solid_linear_thermal_expansion_coefficient)
                      .eval()
                : (crv_u.C * cpv.solid_linear_thermal_expansion_coefficient)
                      .eval();

        local_Jac
            .template block<displacement_size, temperature_size>(
                displacement_index, temperature_index)
            .noalias() -= B.transpose() * uT_coeff * N * w;

        //
        // displacement equation, pressure part (K_up)
        //
        double const up_coeff = cpv.modified_biot;
        Kup.noalias() +=
            B.transpose() * Invariants::identity2 * N * (up_coeff * w);

        //
        // pressure equation, pressure part (K_pp and M_pp).
        //
        double const fluid_density = _ip_data_output[ip].fluid_density;
        double const scaling_factor =
            _process_data.is_volume_balance_equation_type ? 1.0 : fluid_density;

        double const k_rel = fpv.k_rel;
        auto const& K_over_mu = fpv.K_over_mu;
        laplace_p.noalias() +=
            dNdx.transpose() * K_over_mu * dNdx * (k_rel * scaling_factor * w);
        double const dk_rel_dT = crv_T.dk_rel_dT;
        local_Jac
            .template block<pressure_size, temperature_size>(pressure_index,
                                                             temperature_index)
            .noalias() += dNdx.transpose() * K_over_mu * fpv.grad_p * N *
                          (dk_rel_dT * scaling_factor * w);

        double const storage_p_coeff = massBalancePressureStorageCoefficient(
            cpv.vars.porosity, crv_p.fluid_compressibility, cpv.alpha_biot,
            crv_p.beta_SR, has_frozen_liquid_phase ? crv_p.storage_p_fr : 0.0);

        storage_p.noalias() +=
            N.transpose() * N * (storage_p_coeff * scaling_factor * w);

        if (has_frozen_liquid_phase)
        {
            local_Jac
                .template block<pressure_size, temperature_size>(
                    pressure_index, temperature_index)
                .noalias() +=
                N.transpose() * crv_p.J_pT_fr * N * scaling_factor * w;
        }

        laplace_T.noalias() += dNdx.transpose() * fpv.K_pT_thermal_osmosis *
                               dNdx * scaling_factor * w;
        //
        //  RHS, pressure part
        //
        auto const& eps_v_dot = cpv.eps_v_dot;
        double const up_eps_v_dot = up_coeff * eps_v_dot;
        local_rhs.template segment<pressure_size>(pressure_index).noalias() -=
            N * (up_eps_v_dot * scaling_factor * w);

        local_rhs.template segment<pressure_size>(pressure_index).noalias() +=
            dNdx.transpose() * fpv.buoyancy_flux * (scaling_factor * w);
        local_Jac
            .template block<pressure_size, temperature_size>(pressure_index,
                                                             temperature_index)
            .noalias() -= dNdx.transpose() * K_over_mu * b * N *
                          (fluid_density * dk_rel_dT * scaling_factor * w);

        //
        // pressure equation, temperature part (M_pT)
        //

        double const storage_T_coeff = massBalanceTemperatureStorageCoefficient(
            crv_p.beta, has_frozen_liquid_phase ? crv_p.storage_T_fr : 0.0);

        storage_T.noalias() +=
            N.transpose() * storage_T_coeff * N * scaling_factor * w;

        //
        // pressure equation, displacement part.
        //
        // reusing Kup.transpose() if the equation balance type is not volume.

        double const drho_LR_dT = cpv.drho_LR_dT;
        if (!_process_data.is_volume_balance_equation_type)
        {
            Kpu.noalias() += N.transpose() * Invariants::identity2.transpose() *
                             B * (up_coeff * scaling_factor * w);

            //
            // The contribution to Jacobian from d()/ drho
            // drho/dp, d()/ drho drho/dT:
            //
            double const storage_rate = massBalanceStorageRate(
                storage_p_coeff, cpv.p_dot, storage_T_coeff, cpv.T_dot,
                up_eps_v_dot);

            auto const dr_p_drho = massBalanceJacobianWrtFluidDensity(
                ip_data, storage_rate, velocity, fpv.buoyancy_flux);
            local_Jac
                .template block<pressure_size, pressure_size>(pressure_index,
                                                              pressure_index)
                .noalias() += dr_p_drho * crv_p.drho_L_dp;
            local_Jac
                .template block<pressure_size, temperature_size>(
                    pressure_index, temperature_index)
                .noalias() += dr_p_drho * drho_LR_dT;
        }

        //
        // temperature equation
        //
        addEnergyBalanceIntegrationPointTerms(
            energy_balance, ip_data, _ip_data_output[ip], cpv, fpv, crv_T, T);

        // linearized darcy
        double const c_f = crv_T.c_f;
        dKTT_dp_T.noalias() -= N.transpose() * (dNdx * T).transpose() *
                               K_over_mu * dNdx *
                               (fluid_density * c_f * k_rel * w);

        /* TODO (Joerg) Temperature changes due to thermal dilatation of the
         * fluid, which are usually discarded as being very small.
         * Zhou et al. (10.1016/S0020-7683(98)00089-4) states that:
         * "Biot (1956) neglected this term and it is included here for
         * completeness"
         * Keeping the code here in the case these are needed for the named
         * effects in the future.
        if (fluid_compressibility != 0)
        {
            auto const C_el = ip_data.computeElasticTangentStiffness(
                t, x_position, dt, static_cast<double>(T_int_pt));
            auto const solid_skeleton_compressibility =
                1 / solid_material.getBulkModulus(t, x_position, &C_el);
            double const fluid_volumetric_thermal_expansion_coefficient =
                    (fluid_density == 0.0) ? 0.0 : -cpv.drho_LR_dT /
                    fluid_density;
            KTT.noalias() += dNdx.transpose() *
                (-T_int_pt * fluid_volumetric_thermal_expansion_coefficient *
                K_pT_thermal_osmosis / fluid_compressibility) *
                dNdx * w;

            local_rhs.template segment<temperature_size>(temperature_index)
                .noalias() +=
                dNdx.transpose() *
                (-T_int_pt * fluid_volumetric_thermal_expansion_coefficient /
                fluid_compressibility) *
                fluid_density * fpv.k_rel * K_over_mu * b * w;
            MTu part for rhs and Jacobian:
                (-T_int_pt *
                Invariants::trace(solid_linear_thermal_expansion_coefficient) /
                solid_skeleton_compressibility) *
                N.transpose() * identity2.transpose() * B * w;
            KTp part for rhs and Jacobian:
                dNdx.transpose() *
                (T_int_pt * fluid_volumetric_thermal_expansion_coefficient *
                fpv.k_rel * K_over_mu / fluid_compressibility) *
                dNdx * w;
        }
         */
    }

    // temperature equation, pressure part. Before the temperature part below,
    // which consumes the terms.
    local_Jac
        .template block<temperature_size, pressure_size>(temperature_index,
                                                         pressure_index)
        .noalias() += energy_balance.KTp + dKTT_dp_T;

    // temperature equation, temperature part and f_T
    addEnergyBalanceElementTerms(
        std::move(energy_balance),
        local_Jac.template block<temperature_size, temperature_size>(
            temperature_index, temperature_index),
        local_rhs.template segment<temperature_size>(temperature_index), T,
        T_prev, p, dt);

    // displacement equation, pressure part
    local_Jac
        .template block<displacement_size, pressure_size>(displacement_index,
                                                          pressure_index)
        .noalias() -= Kup;

    // pressure equation, temperature part.
    local_Jac
        .template block<pressure_size, temperature_size>(pressure_index,
                                                         temperature_index)
        .noalias() += -storage_T / dt + laplace_T;

    // pressure equation, pressure part.
    local_Jac
        .template block<pressure_size, pressure_size>(pressure_index,
                                                      pressure_index)
        .noalias() += laplace_p + storage_p / dt;

    // pressure equation, displacement part.
    if (_process_data.is_volume_balance_equation_type)
    {
        local_Jac
            .template block<pressure_size, displacement_size>(
                pressure_index, displacement_index)
            .noalias() += Kup.transpose() / dt;
    }
    else
    {
        local_Jac
            .template block<pressure_size, displacement_size>(
                pressure_index, displacement_index)
            .noalias() += Kpu / dt;
    }

    // pressure equation (f_p)
    local_rhs.template segment<pressure_size>(pressure_index).noalias() -=
        laplace_p * p + laplace_T * T + storage_p * (p - p_prev) / dt -
        storage_T * (T - T_prev) / dt;

    // displacement equation (f_u)
    local_rhs.template segment<displacement_size>(displacement_index)
        .noalias() += Kup * p;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
template <typename JacobianBlock, typename ResidualSegment, typename BMatrix,
          typename EffectiveStress>
void ThermoHydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    addMomentumBalanceStiffnessAndResidual(
        JacobianBlock&& jacobian, ResidualSegment&& residual,
        IpData const& ip_data,
        MechanicalPropertyValues<DisplacementDim> const& crv_u,
        BMatrix const& B, EffectiveStress const& effective_stress,
        bool const has_frozen_liquid_phase) const
{
    auto const& N_u = ip_data.N_u;
    double const w = ip_data.integration_weight;
    auto const& b = _process_data.specific_body_force;

    if (has_frozen_liquid_phase)
    {
        residual.noalias() -= B.transpose() * crv_u.r_u_fr * w;
        jacobian.noalias() += B.transpose() * (crv_u.C + crv_u.J_uu_fr) * B * w;
    }
    else
    {
        jacobian.noalias() += B.transpose() * crv_u.C * B * w;
    }

    residual.noalias() -= (B.transpose() * effective_stress -
                           N_u_op(N_u).transpose() * crv_u.rho * b) *
                          w;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
template <typename TemperatureVector>
void ThermoHydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    addEnergyBalanceIntegrationPointTerms(
        EnergyBalanceTerms& energy_balance, IpData const& ip_data,
        IntegrationPointDataForOutput<DisplacementDim> const& ip_data_output,
        CommonPropertyValues<DisplacementDim> const& cpv,
        FlowPropertyValues<DisplacementDim> const& fpv,
        ThermalPropertyValues<DisplacementDim> const& crv_T,
        TemperatureVector const& T) const
{
    auto const& w = ip_data.integration_weight;
    auto const& N = ip_data.N;
    auto const& dNdx = ip_data.dNdx;

    double const T_int_pt = cpv.vars.temperature;
    auto const& velocity = ip_data_output.velocity;
    double const fluid_density = ip_data_output.fluid_density;
    double const c_f = crv_T.c_f;

    //
    // temperature equation, temperature part.
    //
    energy_balance.KTT.noalias() +=
        dNdx.transpose() * crv_T.effective_thermal_conductivity * dNdx * w;
    energy_balance.dKTT_dT_T.noalias() +=
        dNdx.transpose() * crv_T.dlambda_eff_dT * dNdx * T * N * w;

    energy_balance.ip_flux_vector.emplace_back(velocity * fluid_density * c_f);
    // Without any flux correction the flux derivative is as follows. The
    // contribution to KTT is different if any stabilisation scheme is used,
    // but this is ignored for the moment.
    GlobalDimVectorType const dip_flux_vector_dT =
        crv_T.dvelocity_dT * fluid_density * c_f +
        velocity * cpv.drho_LR_dT * c_f;
    energy_balance.dKTT_dT_T.noalias() +=
        N.transpose() * dip_flux_vector_dT.transpose() * dNdx * T * N * w;
    energy_balance.average_velocity_norm += velocity.norm();

    energy_balance.MTT.noalias() += N.transpose() *
                                    (crv_T.sensible_volumetric_heat_capacity -
                                     crv_T.latent_volumetric_heat_capacity) *
                                    N * w;
    energy_balance.dMTT_dT_T.noalias() += N.transpose() * crv_T.J_TT * N * w;

    //
    // temperature equation, pressure part
    //
    energy_balance.KTp.noalias() +=
        dNdx.transpose() * fpv.K_pT_thermal_osmosis * dNdx * (T_int_pt * w);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
template <typename JacobianBlock, typename ResidualSegment,
          typename TemperatureVector, typename TemperaturePrevVector,
          typename PressureVector>
void ThermoHydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    addEnergyBalanceElementTerms(EnergyBalanceTerms&& energy_balance,
                                 JacobianBlock&& jacobian_TT,
                                 ResidualSegment&& residual_T,
                                 TemperatureVector const& T,
                                 TemperaturePrevVector const& T_prev,
                                 PressureVector const& p,
                                 double const dt) const
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    NumLib::assembleAdvectionMatrix(
        _process_data.stabilizer, _ip_data, energy_balance.ip_flux_vector,
        energy_balance.average_velocity_norm /
            static_cast<double>(n_integration_points),
        energy_balance.KTT);

    // temperature equation, temperature part
    jacobian_TT.noalias() += energy_balance.KTT + energy_balance.dKTT_dT_T +
                             energy_balance.MTT / dt + energy_balance.dMTT_dT_T;

    // temperature equation (f_T)
    residual_T.noalias() -= energy_balance.KTT * T +
                            energy_balance.MTT * (T - T_prev) / dt +
                            energy_balance.KTp * p;
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
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::postTimestepConcrete(Eigen::VectorXd const& local_x,
                                           Eigen::VectorXd const& local_x_prev,
                                           double const t, double const dt,
                                           int const /*process_id*/)
{
    auto const T =
        local_x.template segment<temperature_size>(temperature_index);
    // Stored sensible heat integrated over the element, for output.
    double sensible_heat = 0.0;

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto& ip_data = _ip_data[ip];

        auto const space_time_data = makeSpaceTimeData(ip_data.N_u, t, dt);
        IpState const ip_state{ip_data, _ip_data_output[ip]};

        // Evaluated once, before pushBackState(): a later call would add the
        // strain increment of this time step to eps_m a second time.
        auto const cpv = updateCommonProperties(local_x, local_x_prev,
                                                space_time_data, ip_state);

        updateMechanicalProperties(cpv, space_time_data, ip_state);

        ip_data.eps0 =
            ip_data.eps0_prev + (1 - ip_data.phi_fr_prev / cpv.vars.porosity) *
                                    (cpv.eps_prev - ip_data.eps0_prev);
        ip_data.pushBackState();

        // Updates the velocity and viscosity for output.
        updateFlowProperties(cpv, local_x, space_time_data, ip_state);
        auto const hcv = computeHeatCapacities(cpv, space_time_data, ip_state);

        // Stored sensible heat (latent already excluded).
        sensible_heat += hcv.sensible_volumetric_heat_capacity *
                         ip_data.N.dot(T) * ip_data.integration_weight;
    }

    (*_process_data.cell_sensible_heat)[_element.getID()] = sensible_heat;
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

    double phi_fr_avg = 0;
    double fluid_density_avg = 0;
    double viscosity_avg = 0;

    using KV = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    KV sigma_avg = KV::Zero();
    KV sigma_ice_avg = KV::Zero();

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto& ip_data = _ip_data[ip];

        phi_fr_avg += ip_data.phi_fr;
        fluid_density_avg += _ip_data_output[ip].fluid_density;
        viscosity_avg += _ip_data_output[ip].viscosity;
        sigma_avg += ip_data.sigma_eff;
        sigma_ice_avg += ip_data.sigma_eff_ice;
    }

    phi_fr_avg /= n_integration_points;
    fluid_density_avg /= n_integration_points;
    viscosity_avg /= n_integration_points;
    sigma_avg /= n_integration_points;
    sigma_ice_avg /= n_integration_points;

    (*_process_data.element_phi_fr)[_element.getID()] = phi_fr_avg;
    (*_process_data.element_fluid_density)[_element.getID()] =
        fluid_density_avg;
    (*_process_data.element_viscosity)[_element.getID()] = viscosity_avg;

    Eigen::Map<KV>(&(*_process_data.element_stresses)[_element.getID() *
                                                      KV::RowsAtCompileTime]) =
        MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_avg);

    Eigen::Map<KV>(&(
        *_process_data
             .element_ice_stresses)[_element.getID() * KV::RowsAtCompileTime]) =
        MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_ice_avg);

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
