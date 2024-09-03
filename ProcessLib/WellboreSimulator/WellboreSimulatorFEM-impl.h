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

#include <Eigen/Dense>
#include <vector>

#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
namespace WellboreSimulator
{
template <typename ShapeFunction, int GlobalDim>
void WellboreSimulatorFEM<ShapeFunction, GlobalDim>::assemble(
    double const t, double const dt, std::vector<double> const& local_x,
    std::vector<double> const& local_x_prev, std::vector<double>& local_M_data,
    std::vector<double>& local_K_data, std::vector<double>& local_b_data)
{
    auto const local_matrix_size = local_x.size();

    assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

    auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);
    auto local_b = MathLib::createZeroedVector<LocalVectorType>(
        local_b_data, local_matrix_size);

    // Get block matrices
    auto Mvv = local_M.template block<velocity_size, velocity_size>(
        velocity_index, velocity_index);

    auto Mhp = local_M.template block<enthalpy_size, pressure_size>(
        enthalpy_index, pressure_index);
    auto Mhh = local_M.template block<enthalpy_size, enthalpy_size>(
        enthalpy_index, enthalpy_index);

    auto Kpv = local_K.template block<pressure_size, velocity_size>(
        pressure_index, velocity_index);

    auto Kvp = local_K.template block<velocity_size, pressure_size>(
        velocity_index, pressure_index);
    auto Kvv = local_K.template block<velocity_size, velocity_size>(
        velocity_index, velocity_index);

    auto Khh = local_K.template block<enthalpy_size, enthalpy_size>(
        enthalpy_index, enthalpy_index);

    auto Bp = local_b.template segment<pressure_size>(pressure_index);
    auto Bv = local_b.template segment<velocity_size>(velocity_index);
    auto Bh = local_b.template segment<enthalpy_size>(enthalpy_index);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());

    auto const& b = _process_data.specific_body_force;

    MaterialPropertyLib::VariableArray vars;

    // Get material properties
    auto const& medium = *_process_data.media_map.getMedium(_element.getID());
    auto const& liquid_phase = medium.phase("AqueousLiquid");
    auto const& gas_phase = medium.phase("Gas");

    // Get wellbore parameters
    // casing thickness
    auto const t_ca = _process_data.wellbore.casing_thickness(t, pos)[0];
    // wellbore radius
    auto const r_w = _process_data.wellbore.diameter(t, pos)[0] / 2;

    // pipe thickness
    auto const t_p = _process_data.wellbore.pipe_thickness(t, pos)[0];

    // roughness of the wellbore
    auto const xi = _process_data.wellbore.roughness(t, pos)[0];
    // pipe outer radius
    auto const r_o = r_w - t_ca;
    // pipe inner radius
    auto const r_i = r_o - t_p;

    // get reservoir properties
    NodalVectorType T_r =
        _process_data.reservoir_properties.temperature.getNodalValuesOnElement(
            _element, t);
    NodalVectorType p_r =
        _process_data.reservoir_properties.pressure.getNodalValuesOnElement(
            _element, t);
    NodalVectorType PI =
        _process_data.productivity_index.getNodalValuesOnElement(_element, t);
    auto const k_r =
        _process_data.reservoir_properties.thermal_conductivity(t, pos)[0];
    auto const rho_r = _process_data.reservoir_properties.density(t, pos)[0];
    auto const c_r =
        _process_data.reservoir_properties.specific_heat_capacity(t, pos)[0];

    for (unsigned ip(0); ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);

        auto& ip_data = _ip_data[ip];
        auto const& N = ip_data.N;
        auto const& dNdx = ip_data.dNdx;
        auto const& w = ip_data.integration_weight;
        auto& mix_density = ip_data.mix_density;
        auto& temperature = ip_data.temperature;
        auto& steam_mass_frac = ip_data.dryness;
        auto& vapor_volume_frac = ip_data.vapor_volume_fraction;
        auto& vapor_mass_flowrate = ip_data.vapor_mass_flow_rate;
        auto& liquid_mass_flowrate = ip_data.liquid_mass_flow_rate;

        double p_int_pt = 0.0;
        double v_int_pt = 0.0;
        double h_int_pt = 0.0;

        NumLib::shapeFunctionInterpolate(local_x, N, p_int_pt, v_int_pt,
                                         h_int_pt);

        double p_prev_int_pt = 0.0;
        double v_prev_int_pt = 0.0;
        double h_prev_int_pt = 0.0;

        NumLib::shapeFunctionInterpolate(local_x_prev, N, p_prev_int_pt,
                                         v_prev_int_pt, h_prev_int_pt);

        double vdot_int_pt = (v_int_pt - v_prev_int_pt) / dt;

        // calculate fluid properties

        const double pi = std::numbers::pi;

        vars.liquid_phase_pressure = p_int_pt;
        vars.enthalpy = h_int_pt;

        double liquid_water_density =
            liquid_phase
                .property(MaterialPropertyLib::PropertyType::saturation_density)
                .template value<double>(vars, pos, t, dt);
        double const vapour_water_density =
            gas_phase
                .property(MaterialPropertyLib::PropertyType::saturation_density)
                .template value<double>(vars, pos, t, dt);

        double const h_sat_liq_w =
            liquid_phase
                .property(
                    MaterialPropertyLib::PropertyType::saturation_enthalpy)
                .template value<double>(vars, pos, t, dt);
        double const h_sat_vap_w =
            gas_phase
                .property(
                    MaterialPropertyLib::PropertyType::saturation_enthalpy)
                .template value<double>(vars, pos, t, dt);

        // TODO add a function to calculate dryness with constrain of 0
        // to 1.
        double const dryness = std::max(
            0., (h_int_pt - h_sat_liq_w) / (h_sat_vap_w - h_sat_liq_w));
        steam_mass_frac = dryness;

        double const T_int_pt =
            (dryness == 0)
                ? liquid_phase
                      .property(MaterialPropertyLib::PropertyType::temperature)
                      .template value<double>(vars, pos, t, dt)
                : gas_phase
                      .property(MaterialPropertyLib::PropertyType::
                                    saturation_temperature)
                      .template value<double>(vars, pos, t, dt);
        temperature = T_int_pt;
        vars.temperature = T_int_pt;

        // For the calculation of the void fraction of vapour,
        // see Rohuani, Z., and E. Axelsson. "Calculation of volume void
        // fraction in a subcooled and quality region." International
        // Journal of Heat and Mass Transfer 17 (1970): 383-393.

        // profile parameter of drift flux
        double C_0 = 1 + 0.12 * (1 - dryness);

        // For the surface tension calculation, see
        // Cooper, J. R., and R. B. Dooley. "IAPWS release on surface
        // tension of ordinary water substance." International Association
        // for the Properties of Water and Steam (1994).
        double const sigma_gl = 0.2358 *
                                std::pow((1 - T_int_pt / 647.096), 1.256) *
                                (1 - 0.625 * (1 - T_int_pt / 647.096));
        // drift flux velocity
        double const u_gu =
            1.18 * (1 - dryness) *
            std::pow((9.81) * sigma_gl *
                         (liquid_water_density - vapour_water_density),
                     0.25) /
            std::pow(liquid_water_density, 0.5);

        // solving void fraction of vapor: Rouhani-Axelsson
        double alpha = 0;
        if (dryness != 0)
        {
            // Local Newton solver
            using LocalJacobianMatrix =
                Eigen::Matrix<double, 1, 1, Eigen::RowMajor>;
            using LocalResidualVector = Eigen::Matrix<double, 1, 1>;
            using LocalUnknownVector = Eigen::Matrix<double, 1, 1>;
            LocalJacobianMatrix J_loc;

            Eigen::PartialPivLU<LocalJacobianMatrix> linear_solver(1);

            auto const update_residual = [&](LocalResidualVector& residual)
            {
                calculateResidual(alpha, vapour_water_density,
                                  liquid_water_density, v_int_pt, dryness, C_0,
                                  u_gu, residual);
            };

            auto const update_jacobian = [&](LocalJacobianMatrix& jacobian)
            {
                calculateJacobian(
                    alpha, vapour_water_density, liquid_water_density, v_int_pt,
                    dryness, C_0, u_gu,
                    jacobian);  // for solution dependent Jacobians
            };

            auto const update_solution =
                [&](LocalUnknownVector const& increment)
            {
                // increment solution vectors
                alpha += increment[0];
            };

            const int maximum_iterations(20);
            const double residuum_tolerance(1.e-10);
            const double increment_tolerance(0);

            auto newton_solver = NumLib::NewtonRaphson(
                linear_solver, update_jacobian, update_residual,
                update_solution,
                {maximum_iterations, residuum_tolerance, increment_tolerance});

            auto const success_iterations = newton_solver.solve(J_loc);

            if (!success_iterations)
            {
                WARN(
                    "Attention! Steam void fraction has not been correctly "
                    "calculated!");
            }
        }

        vapor_volume_frac = alpha;

        if (alpha == 0)
        {
            liquid_water_density =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::density)
                    .template value<double>(vars, pos, t, dt);
        }

        mix_density =
            vapour_water_density * alpha + liquid_water_density * (1 - alpha);

        auto& mix_density_prev = ip_data.mix_density_prev;
        vars.density = mix_density;

        auto const rho_dot = (mix_density - mix_density_prev) / dt;

        double const liquid_water_velocity_act =
            (alpha == 0) ? v_int_pt
                         : (1 - dryness) * mix_density * v_int_pt /
                               (1 - alpha) / liquid_water_density;
        double const vapor_water_velocity_act =
            (alpha == 0) ? 0
                         : dryness * mix_density * v_int_pt /
                               (alpha * vapour_water_density);

        vapor_mass_flowrate = vapor_water_velocity_act * vapour_water_density *
                              pi * r_i * r_i * alpha;

        liquid_mass_flowrate = liquid_water_velocity_act *
                               liquid_water_density * pi * r_i * r_i *
                               (1 - alpha);

        // Slip parameter between two phases,
        // see Akbar, Somaieh, N. Fathianpour, and Rafid Al Khoury. "A finite
        // element model for high enthalpy two-phase flow in geothermal
        // wellbores." Renewable Energy 94 (2016): 223-236.
        double const gamma =
            alpha * liquid_water_density * vapour_water_density * mix_density /
            (1 - alpha) /
            std::pow((alpha * C_0 * vapour_water_density +
                      (1 - alpha * C_0) * liquid_water_density),
                     2) *
            std::pow((C_0 - 1) * v_int_pt + u_gu, 2);

        double const miu =
            liquid_phase.property(MaterialPropertyLib::PropertyType::viscosity)
                .template value<double>(vars, pos, t, dt);
        double const Re = mix_density * v_int_pt * 2 * r_i / miu;

        // Wall friction coefficient,
        // Musivand Arzanfudi, Mehdi, and Rafid Al‐Khoury. "A compressible
        // two‐fluid multiphase model for CO2 leakage through a wellbore."
        // International Journal for Numerical Methods in Fluids 77.8 (2015):
        // 477-507.
        double f = 0.0;
        if (Re > 10 && Re <= 2400)
            f = 16 / Re;
        else if (Re > 2400)
            f = std::pow(std::log(xi / 3.7 / r_i) -
                             5.02 / Re * std::log(xi / 3.7 / r_i + 13 / Re),
                         -2) /
                16;

        double Q_hx = 0;
        double const T_r_int_pt = N.dot(T_r);
        // conductive heat exchange between wellbore and formation
        if (_process_data.has_heat_exchange_with_formation)
        {
            // See Zhang, Pan, Pruess, Finsterle (2011). A time-convolution
            // approach for modeling heat exchange between a wellbore and
            // surrounding formation. Geothermics 40, 261-266.
            const double alpha_r = k_r / rho_r / c_r;
            const double t_d = alpha_r * t / (r_i * r_i);

            double beta;
            if (t_d < 2.8)
                beta = std::pow((pi * t_d), -0.5) + 0.5 -
                       0.25 * std::pow((t_d / pi), 0.5) + 0.125 * t_d;
            else
                beta = 2 * (1 / (std::log(4 * t_d) - 2 * 0.57722) -
                            0.57722 /
                                std::pow((std::log(4 * t_d) - 2 * 0.57722), 2));

            const double P_c = 2 * pi * r_i;
            Q_hx = P_c * k_r * (T_r_int_pt - T_int_pt) / r_i * beta;
        }

        // mass exchange with reservoir
        double const p_r_int_pt = N.dot(p_r);
        double const PI_int_pt = N.dot(PI);
        double Q_mx = PI_int_pt * (p_int_pt - p_r_int_pt);

        // advective momentum and energy exchange with reservoir due to the mass
        // exchange
        double Q_mom = 0;
        double Q_ene = 0;
        if (Q_mx != 0)
        {
            Q_mom = Q_mx * v_int_pt;
            // Only single-phase liquid condition is considered now
            // TODO: update the two-phase flowing enthalpy from the feed zone.
            vars.liquid_phase_pressure = p_r_int_pt;
            vars.temperature = T_r_int_pt;
            double const h_fres =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::enthalpy)
                    .template value<double>(vars, pos, t, dt);
            Q_ene = Q_mx * h_fres;
        }

        // M matrix assembly
        Mvv.noalias() += w * N.transpose() * mix_density * N;

        Mhp.noalias() += -w * N.transpose() * N;
        Mhh.noalias() += w * N.transpose() * mix_density * N;

        // K matrix assembly
        Kpv.noalias() += w * dNdx.transpose() * N * mix_density;

        Kvp.noalias() += w * N.transpose() * dNdx;
        Kvv.noalias() += w * N.transpose() * rho_dot * N;

        Khh.noalias() += w * N.transpose() * mix_density * v_int_pt * dNdx;

        // b matrix assembly
        Bp.noalias() += w * N.transpose() * rho_dot + w * N.transpose() * Q_mx;

        Bv.noalias() +=
            w * dNdx.transpose() * mix_density * v_int_pt * v_int_pt +
            w * dNdx.transpose() * gamma -
            w * N.transpose() * f * mix_density * std::abs(v_int_pt) *
                v_int_pt / (4 * r_i) -
            w * N.transpose() * Q_mom;

        Bh.noalias() +=
            -1 / 2 * w * N.transpose() * rho_dot * v_int_pt * v_int_pt -
            w * N.transpose() * mix_density * v_int_pt * vdot_int_pt +
            1 / 2 * w * dNdx.transpose() * mix_density * v_int_pt * v_int_pt *
                v_int_pt +
            w * N.transpose() * (Q_hx / pi / r_i / r_i) -
            w * N.transpose() * Q_ene;

        if (_process_data.has_gravity)
        {
            NodalVectorType gravity_operator =
                N.transpose() * b * w * _element_direction[2];

            Bv.noalias() += gravity_operator * mix_density;
            Bh.noalias() += gravity_operator * mix_density * v_int_pt;
        }
    }

    //     debugging
    //    std::string sep = "\n----------------------------------------\n";
    //    Eigen::IOFormat CleanFmt(6, 0, ", ", "\n", "[", "]");
    //    std::cout << local_M.format(CleanFmt) << sep;
    //    std::cout << local_K.format(CleanFmt) << sep;
    //    std::cout << local_b.format(CleanFmt) << sep;
}

}  // namespace WellboreSimulator
}  // namespace ProcessLib
