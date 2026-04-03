// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "SaturationLuMcCartney.h"

#include <algorithm>
#include <cmath>
#include <numbers>

#include "MaterialLib/PhysicalConstant.h"
using namespace boost::math::differentiation;

#include "MaterialLib/MPL/Medium.h"

namespace MaterialPropertyLib
{

template <typename PressureType, typename TemperatureType>
promote<PressureType, TemperatureType> SaturationLuMcCartney::water_content(
    PressureType p_cap, TemperatureType T,
    VariableArray const& variable_array) const
{
    return adsorptive_water_content(p_cap, T, variable_array) +
           capillary_water_content(p_cap, T, variable_array);
}

template <typename PressureType, typename TemperatureType>
promote<PressureType, TemperatureType>
SaturationLuMcCartney::adsorptive_water_content(
    PressureType p_cap, TemperatureType T,
    VariableArray const& variable_array) const
{
    return theta_a_max(T, variable_array) *
           (1 - pow(exp((p_cap - psi_max(T)) / p_cap), M));
}

template <typename PressureType, typename TemperatureType>
promote<PressureType, TemperatureType>
SaturationLuMcCartney::capillary_water_content(
    PressureType p_cap, TemperatureType T,
    VariableArray const& variable_array) const
{
    auto const theta_mean =
        (variable_array.porosity -
         adsorptive_water_content(p_cap, T, variable_array)) /
        2.0;
    auto const chi_frac = (chi(T) + T) / (chi_r + Tr);
    auto const erf_term =
        erf(sqrt(2) * (p_cap / psi_c(T, variable_array) * chi_frac - 1));
    auto const power_term = pow(alpha(T) * p_cap * chi_frac, N);
    auto const cap_term =
        theta_mean * (1 - erf_term) * pow(1 + power_term, ((1 / N) - 1));
    return cap_term;
}

template <typename TemperatureType>
promote<TemperatureType> SaturationLuMcCartney::theta_a_max(
    TemperatureType T, VariableArray const& variable_array) const
{
    return (1 - variable_array.porosity) * (CEC(T) / zeta_s(T) + bw);
}

template <typename TemperatureType>
promote<TemperatureType> SaturationLuMcCartney::psi_max(TemperatureType T) const
{
    return MaterialLib::PhysicalConstant::IdealGasConstant * T * c(T) /
           (3.0 * nu_w);
}

template <typename TemperatureType>
promote<TemperatureType> SaturationLuMcCartney::c(TemperatureType T) const
{
    return exp((E1_minus_EL) /
               (MaterialLib::PhysicalConstant::IdealGasConstant * T));
}

template <typename TemperatureType>
promote<TemperatureType> SaturationLuMcCartney::CEC(TemperatureType T) const
{
    auto const cos_arg = A * std::numbers::pi * (B * T - T1) / (T2 - T1);
    return 0.5 * CEC_max * (cos(cos_arg) + 1);
}

template <typename TemperatureType>
promote<TemperatureType> SaturationLuMcCartney::psi_c(
    TemperatureType T, VariableArray const& variable_array) const
{
    return A_H(T) / (6 * std::numbers::pi) *
           pow(theta_a_max(T, variable_array) / (SSA * variable_array.density),
               -3.0);
}

template <typename TemperatureType>
promote<TemperatureType> SaturationLuMcCartney::A_H(TemperatureType T) const
{
    auto const dielectric_term =
        pow((epsilon_s - epsilon_w(T)) / (epsilon_s + epsilon_w(T)), 2);
    auto const RI_term =
        pow(n_s * n_s - n_w * n_w, 2) / pow(n_s * n_s + n_w * n_w, 3. / 2.);
    auto const prefactor_dielectic =
        (3 * MaterialLib::PhysicalConstant::BoltzmannConstant * T) / 4.;
    auto const prefactor_RI =
        (3 * MaterialLib::PhysicalConstant::PlanckConstant * nu_e) /
        (16 * sqrt(2));
    return prefactor_dielectic * dielectric_term + prefactor_RI * RI_term;
}

template <typename TemperatureType>
promote<TemperatureType> SaturationLuMcCartney::chi(TemperatureType T) const
{
    return delta_h(T) / C1;
}

template <typename TemperatureType>
promote<TemperatureType> SaturationLuMcCartney::delta_h(TemperatureType T) const
{
    // delta_hr = -0.516; //J / m ^ 2 high plasticity clay
    auto const delta_hr = -0.5151936137548038;
    return delta_hr * pow((1 - Tr) / (1 - T), 0.38);
}

template <typename TemperatureType>
promote<TemperatureType> SaturationLuMcCartney::alpha(TemperatureType T) const
{
    auto const psi_aev = alpha_0_inv * exp(eta_alpha * (Tr / T - 1));
    return 1.0 / psi_aev;
}

template <typename TemperatureType>
promote<TemperatureType> SaturationLuMcCartney::epsilon_w(
    TemperatureType T) const
{
    auto const r_w = rho_w(T) / 1000.0;
    return pow(
        10, (0.7017 + 642.0 / T - 1.167e5 / pow(T, 2) + 9.190e6 / pow(T, 3) +
             (1.667 - 11.41 / T - 3.526e4 / pow(T, 2)) * log(r_w) / log(10)) +
                1);
}

template <typename TemperatureType>
promote<TemperatureType> SaturationLuMcCartney::zeta_s(TemperatureType T) const
{
    return (c(T) - 1) / c(T) * CEC(T) / nu_mr;
}

template <typename TemperatureType>
promote<TemperatureType> SaturationLuMcCartney::rho_w(TemperatureType T) const
{
    double const pi = 0.0;
    auto const tau = ref_T_ / T;

    return ref_p_ /
           (MaterialLib::PhysicalConstant::SpecificGasConstant::WaterVapour *
            T * gibbs_free_energy_.get_dgamma_dpi(tau, pi));
}

SaturationLuMcCartney::SaturationLuMcCartney(std::string name,
                                             std::string const& material)
    : material_(material),
      gibbs_free_energy_(
          MaterialLib::Fluid::DimensionLessGibbsFreeEnergyRegion1())
{
    name_ = std::move(name);
    if (material_ == "MX80")
    {
        E1_minus_EL = 7533.930402019812;
        nu_mr = 2.752035e-01;
        CEC_max = 0.8004082526701221;
        bw = 0.18;
        M = 0.105;
        N = 1.15;
        SSA = 700.0e3;
        epsilon_s = 3.2089251122574;
        n_s = 1.3298392371979677;
        eta_alpha = 5.0;
        alpha_0_inv = 3.3e6;
        chi_r = delta_h(Tr) / C1;
    }
    else if (material_ == "BoomClay")
    {
        E1_minus_EL = 6803;
        nu_mr = 7.994125e-02;
        CEC_max = 0.2534798545419938;
        bw = 0.11;
        M = 0.085;
        N = 1.29;
        SSA = 260e3;
        epsilon_s = 1.3267932723717908;
        n_s = 35.435652398543425;
        eta_alpha = 0.6;
        alpha_0_inv = 36e3;
        chi_r = delta_h(Tr) / C1;
    }
    else if (material_ == "FEBEX")
    {
        E1_minus_EL = 7904.547103777231;
        nu_mr = 3.374814e-01;
        CEC_max = 1.0108839312584308;
        bw = 0.18;
        M = 0.132;
        N = 1.22;
        SSA = 860.0e3;
        epsilon_s = 0.940367803768519;
        n_s = 748.7804881357005;
        eta_alpha = 3.5;
        alpha_0_inv = 330e3;
        chi_r = delta_h(Tr) / C1;
    }
    else if (material_ == "GMZ01")
    {
        E1_minus_EL = 8369.804359014617;
        nu_mr = 1.505956e-01;
        CEC_max = 0.8004082526701221;
        bw = 0.25;
        M = 0.135;
        N = 1.3;
        SSA = 700.0e3;
        epsilon_s = 2.1590265530454893;
        n_s = 2175.885110079861;
        eta_alpha = 8.5;
        alpha_0_inv = 8.0e6;
        chi_r = delta_h(Tr) / C1;
    }

    else
    {
        OGS_FATAL(
            "Material '{}' not implemented in SaturationLuMcCartney model.",
            material_);
    }
}

PropertyDataType SaturationLuMcCartney::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    const double p_cap = variable_array.capillary_pressure;
    const double T = variable_array.temperature;

    auto const& medium = *std::get<Medium*>(scale_);
    auto const& solid_phase = medium.phase(PhaseName::Solid);
    auto const& porosity_property = medium[PropertyType::porosity];
    auto const& solid_denisty_property = solid_phase[PropertyType::density];

    auto const phi =
        std::get<double>(porosity_property.value(variable_array, pos, t, dt));
    auto const rho_s = std::get<double>(
        solid_denisty_property.value(variable_array, pos, t, dt));
    VariableArray rho_d_and_phi;
    rho_d_and_phi.density = rho_s * (1.0 - phi);
    rho_d_and_phi.porosity = phi;
    double const p_cut = 1e-10;

    auto const S_L_max = water_content(p_cut, T, rho_d_and_phi) / phi;
    auto const p_cap_max = psi_max(T);
    if (p_cap <= p_cut)
    {
        return S_L_max;
    }
    if (p_cap >= p_cap_max)
    {
        return 0.0;
    }
    return water_content(p_cap, T, rho_d_and_phi) / phi;
}

PropertyDataType SaturationLuMcCartney::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    if (variable != Variable::capillary_pressure &&
        variable != Variable::temperature)
    {
        OGS_FATAL(
            "SaturationLuMcCartney::dValue is implemented for derivatives "
            "with respect to capillary pressure or temperature only.");
    }

    const double p_cap = variable_array.capillary_pressure;
    const double T = variable_array.temperature;

    auto const& medium = *std::get<Medium*>(scale_);
    auto const& solid_phase = medium.phase(PhaseName::Solid);
    auto const& porosity_property = medium[PropertyType::porosity];
    auto const& solid_denisty_property = solid_phase[PropertyType::density];

    auto const phi =
        std::get<double>(porosity_property.value(variable_array, pos, t, dt));
    auto const rho_s = std::get<double>(
        solid_denisty_property.value(variable_array, pos, t, dt));
    VariableArray rho_d_and_phi;
    rho_d_and_phi.density = rho_s * (1.0 - phi);
    rho_d_and_phi.porosity = phi;

    auto const p_cap_max = psi_max(T);
    auto const p_cut = 1.0e-10;

    auto const variables = make_ftuple<double, 1, 1>(p_cap, T);
    auto const& v = std::get<0>(variables);
    auto const& w = std::get<1>(variables);
    auto const y = water_content(v, w, rho_d_and_phi);
    double DS = 0.0;
    if (p_cap >= p_cap_max)
    {
        return 0.0;
    }
    if (variable == Variable::capillary_pressure)
    {
        if (p_cap <= 1.0e-10)
        {
            return 0.;
        }
        return y.derivative(1, 0) / phi;
    }
    if (variable == Variable::temperature)
    {
        if (p_cap <= p_cut)
        {
            auto const variables_b2 = make_ftuple<double, 1, 1>(p_cut, T);
            auto const& v_b2 = std::get<0>(variables_b2);
            auto const& w_b2 = std::get<1>(variables_b2);
            auto const y_b2 = water_content(v_b2, w_b2, rho_d_and_phi);
            return y_b2.derivative(0, 1) / phi;
        }
        DS = y.derivative(0, 1);
    }
    return DS / phi;
}

PropertyDataType SaturationLuMcCartney::d2Value(
    VariableArray const& variable_array, Variable const variable1,
    Variable const variable2, ParameterLib::SpatialPosition const& pos,
    double const t, double const dt) const
{
    (void)variable1;
    (void)variable2;
    assert(((variable1 == Variable::capillary_pressure) ||
            (variable1 == Variable::temperature)) &&
           ((variable2 == Variable::capillary_pressure) ||
            (variable2 == Variable::temperature)) &&
           "SaturationLuMcCartney::d2Value is implemented for  derivatives "
           "with respect to capillary pressure only.");

    const double p_cap = variable_array.capillary_pressure;
    const double T = variable_array.temperature;

    auto const& medium = *std::get<Medium*>(scale_);
    auto const& solid_phase = medium.phase(PhaseName::Solid);
    auto const& porosity_property = medium[PropertyType::porosity];
    auto const& solid_denisty_property = solid_phase[PropertyType::density];

    auto const phi =
        std::get<double>(porosity_property.value(variable_array, pos, t, dt));
    auto const rho_s = std::get<double>(
        solid_denisty_property.value(variable_array, pos, t, dt));
    VariableArray rho_d_and_phi;
    rho_d_and_phi.density = rho_s * (1.0 - phi);
    rho_d_and_phi.porosity = phi;

    auto const p_cap_max = psi_max(T);
    auto const p_cut = 1e-10;

    if (p_cap >= p_cap_max)
    {
        return 0.;
    }
    if (p_cap <= p_cut && (variable1 == Variable::capillary_pressure ||
                           variable2 == Variable::capillary_pressure))
    {
        return 0.;
    }

    auto const variables = make_ftuple<double, 2, 2>(p_cap, T);
    auto const& v = std::get<0>(variables);  // Up to Nw derivatives at w=11
    auto const& w = std::get<1>(variables);
    auto const y = water_content(v, w, rho_d_and_phi);
    double DDS = 0.0;
    auto returnVariableTuple = [](const Variable v1,
                                  const Variable v2) -> std::tuple<int, int>
    {
        if (v1 == Variable::capillary_pressure &&
            v2 == Variable::capillary_pressure)
        {
            return {2, 0};
        }
        if (v1 == Variable::temperature && v2 == Variable::temperature)
        {
            return {0, 2};
        }
        return {1, 1};
    };
    auto [i, j] = returnVariableTuple(variable1, variable2);
    if (j > 0)
    {
        if (p_cap <= p_cut)
        {
            auto const variables_b2 = make_ftuple<double, 2, 2>(p_cut, T);
            auto const& v_b2 = std::get<0>(variables_b2);
            auto const& w_b2 = std::get<1>(variables_b2);
            auto const y_b2 = water_content(v_b2, w_b2, rho_d_and_phi);
            return y_b2.derivative(i, j) / phi;
        }
    }
    DDS = y.derivative(i, j);
    return DDS / phi;
}
}  // namespace MaterialPropertyLib
