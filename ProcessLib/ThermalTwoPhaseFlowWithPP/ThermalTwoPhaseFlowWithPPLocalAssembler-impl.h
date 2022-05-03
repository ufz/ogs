/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**
 * Common convenitions for naming:
 * x_air_nonwet           mass fraction of gas component(e.g air) in nonwetting
 * phase (gas component doesn't include water vapour, same for the following)
 * x_air_nonwet           molar fraction of gas component in nonwetting phase
 * x_water_nonwet         molar fraction of vapour in nonwetting phase
 * p_vapour_nonwet         water vapour pressure
 * p_gas_nonwet           partial pressure of gas component
 * mol_density_nonwet         molar density of nonwetting phase
 * mol_density_water          molar density of water
 * density_water              mass density of water
 * density_air_nonwet         mass density of gas component in the nonwetting
 * phase density_nonwet_vapour       mass density of vapour in the nonwetting
 * phase density_nonwet             mass density of the nonwetting phase
 * density_wet                mass density of wetting pahse
 * density_solid              mass density of the solid phase
 * velocity_nonwet              velocity of nonwetting phase
 * velocity_wet                 velocity of wetting phase
 * heat_capacity_dry_gas        heat capacity of dry gas
 * heat_capacity_water_vapour    heat capacity of water vapour
 * heat_capacity_water          heat capacity of liquid water
 * heat_capacity_solid          heat capacity of soil grain
 * latent_heat_evaporation      latent heat for evaporation(water to vapour)
 * enthalpy_nonwet_gas          enthalpy of gas component in the nonwetting
 * phase enthalpy_nonwet_vapour        enthalpy of water vapour in the
 * nonwetting phase enthalpy_wet                 enthalpy of wetting phase
 * enthalpy_nonwet                 enthalpy of the nonwetting phase
 * internal_energy_nonwet        specific internal energy for the nonwetting
 * phase internal_energy_wet           specific internal energy for the wetting
 * phase heat_conductivity_dry_solid   heat conductivity of the dry porous
 * medium heat_conductivity_wet_solid   heat conductivity of the fully saturated
 * porous medium heat_conductivity_unsaturated   heat conductivity of the
 * unsaturated porous medium
 */
#pragma once

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/Utils/FormEffectiveThermalConductivity.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/Function/Interpolation.h"
#include "ThermalTwoPhaseFlowWithPPLocalAssembler.h"
#include "ThermalTwoPhaseFlowWithPPProcessData.h"

namespace ProcessLib
{
namespace ThermalTwoPhaseFlowWithPP
{
template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
void ThermalTwoPhaseFlowWithPPLocalAssembler<
    ShapeFunction, IntegrationMethod,
    GlobalDim>::assemble(double const t, double const dt,
                         std::vector<double> const& local_x,
                         std::vector<double> const& /*local_xdot*/,
                         std::vector<double>& local_M_data,
                         std::vector<double>& local_K_data,
                         std::vector<double>& local_b_data)
{
    using MaterialLib::PhysicalConstant::CelsiusZeroInKelvin;
    using MaterialLib::PhysicalConstant::IdealGasConstant;

    auto const local_matrix_size = local_x.size();

    assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

    auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);
    auto local_b = MathLib::createZeroedVector<LocalVectorType>(
        local_b_data, local_matrix_size);

    auto Map =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Mapc = local_M.template block<nonwet_pressure_size, cap_pressure_size>(
        nonwet_pressure_matrix_index, cap_pressure_matrix_index);
    auto Max =
        local_M.template block<nonwet_pressure_size, tot_mol_frac_contam_size>(
            nonwet_pressure_matrix_index, tot_mol_frac_contam_matrix_index);
    auto Mat = local_M.template block<nonwet_pressure_size, temperature_size>(
        nonwet_pressure_matrix_index, temperature_matrix_index);

    auto Mwp = local_M.template block<cap_pressure_size, nonwet_pressure_size>(
        cap_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Mwpc = local_M.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, cap_pressure_matrix_index);
    auto Mwx =
        local_M.template block<cap_pressure_size, tot_mol_frac_contam_size>(
            cap_pressure_matrix_index, tot_mol_frac_contam_matrix_index);
    auto Mwt = local_M.template block<cap_pressure_size, temperature_size>(
        cap_pressure_matrix_index, temperature_matrix_index);

    auto Mcp =
        local_M.template block<tot_mol_frac_contam_size, nonwet_pressure_size>(
            tot_mol_frac_contam_matrix_index, nonwet_pressure_matrix_index);
    auto Mcpc =
        local_M.template block<tot_mol_frac_contam_size, cap_pressure_size>(
            tot_mol_frac_contam_matrix_index, cap_pressure_matrix_index);
    auto Mcx = local_M.template block<tot_mol_frac_contam_size,
                                      tot_mol_frac_contam_size>(
        tot_mol_frac_contam_matrix_index, tot_mol_frac_contam_matrix_index);
    auto Mct =
        local_M.template block<tot_mol_frac_contam_size, temperature_size>(
            tot_mol_frac_contam_matrix_index, temperature_matrix_index);

    auto Mep = local_M.template block<temperature_size, nonwet_pressure_size>(
        temperature_matrix_index, nonwet_pressure_matrix_index);
    auto Mepc = local_M.template block<temperature_size, cap_pressure_size>(
        temperature_matrix_index, cap_pressure_matrix_index);
    auto Mex =
        local_M.template block<temperature_size, tot_mol_frac_contam_size>(
            temperature_matrix_index, tot_mol_frac_contam_matrix_index);
    auto Met = local_M.template block<temperature_size, temperature_size>(
        temperature_matrix_index, temperature_matrix_index);

    NodalMatrixType laplace_operator =
        NodalMatrixType::Zero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

    auto Kap =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Kapc = local_K.template block<nonwet_pressure_size, cap_pressure_size>(
        nonwet_pressure_matrix_index, cap_pressure_matrix_index);
    auto Kax =
        local_K.template block<nonwet_pressure_size, tot_mol_frac_contam_size>(
            nonwet_pressure_matrix_index, tot_mol_frac_contam_matrix_index);
    auto Kat = local_K.template block<nonwet_pressure_size, temperature_size>(
        nonwet_pressure_matrix_index, temperature_matrix_index);

    auto Kwp = local_K.template block<cap_pressure_size, nonwet_pressure_size>(
        cap_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Kwpc = local_K.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, cap_pressure_matrix_index);
    auto Kwx =
        local_K.template block<cap_pressure_size, tot_mol_frac_contam_size>(
            cap_pressure_matrix_index, tot_mol_frac_contam_matrix_index);
    auto Kwt = local_K.template block<cap_pressure_size, temperature_size>(
        cap_pressure_matrix_index, temperature_matrix_index);

    auto Kcp =
        local_K.template block<tot_mol_frac_contam_size, nonwet_pressure_size>(
            tot_mol_frac_contam_matrix_index, nonwet_pressure_matrix_index);
    auto Kcpc =
        local_K.template block<tot_mol_frac_contam_size, cap_pressure_size>(
            tot_mol_frac_contam_matrix_index, cap_pressure_matrix_index);
    auto Kcx = local_K.template block<tot_mol_frac_contam_size,
                                      tot_mol_frac_contam_size>(
        tot_mol_frac_contam_matrix_index, tot_mol_frac_contam_matrix_index);
    auto Kct =
        local_K.template block<tot_mol_frac_contam_size, temperature_size>(
            tot_mol_frac_contam_matrix_index, temperature_matrix_index);

    auto Kep = local_K.template block<temperature_size, nonwet_pressure_size>(
        temperature_matrix_index, nonwet_pressure_matrix_index);
    auto Kepc = local_K.template block<temperature_size, cap_pressure_size>(
        temperature_matrix_index, cap_pressure_matrix_index);
    auto Ket = local_K.template block<temperature_size, temperature_size>(
        temperature_matrix_index, temperature_matrix_index);

    auto Ba = local_b.template segment<nonwet_pressure_size>(
        nonwet_pressure_matrix_index);
    auto Bw =
        local_b.template segment<cap_pressure_size>(cap_pressure_matrix_index);
    auto Bc = local_b.template segment<tot_mol_frac_contam_size>(
        tot_mol_frac_contam_matrix_index);
    auto Be =
        local_b.template segment<temperature_size>(temperature_matrix_index);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());

    auto const num_nodes = ShapeFunction::NPOINTS;
    auto const pg_nodal_values =
        Eigen::Map<const NodalVectorType>(&local_x[0], num_nodes);
    auto const pc_nodal_values =
        Eigen::Map<const NodalVectorType>(&local_x[num_nodes], num_nodes);

    MaterialPropertyLib::VariableArray vars;

    GlobalDimMatrixType const& I(
        GlobalDimMatrixType::Identity(GlobalDim, GlobalDim));

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& ip_data = _ip_data[ip];
        auto const& N = ip_data.N;
        auto const& dNdx = ip_data.dNdx;
        auto const& w = ip_data.integration_weight;
        auto const& mass_operator = ip_data.mass_operator;
        auto const& diffusion_operator = ip_data.diffusion_operator;
        double pg_int_pt = 0.;
        double pc_int_pt = 0.;
        double Xc_int_pt = 0.;
        double T_int_pt = 0.;
        NumLib::shapeFunctionInterpolate(local_x, N, pg_int_pt, pc_int_pt,
                                         Xc_int_pt, T_int_pt);

        _pressure_wetting[ip] = pg_int_pt - pc_int_pt;
        double const ideal_gas_constant_times_T_int_pt =
            IdealGasConstant * T_int_pt;
        vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
            T_int_pt;
        vars[static_cast<int>(
            MaterialPropertyLib::Variable::capillary_pressure)] = pc_int_pt;
        vars[static_cast<int>(MaterialPropertyLib::Variable::phase_pressure)] =
            pg_int_pt;

        auto const& medium =
            *_process_data.media_map->getMedium(this->_element.getID());
        auto const& liquid_phase = medium.phase("AqueousLiquid");
        auto const& solid_phase = medium.phase("Solid");
        auto const& gas_phase = medium.phase("Gas");

        auto const& water_vapour_component = gas_phase.component("w");
        auto const& dry_air_component = gas_phase.component("a");
        auto const& contam_vapour_component = gas_phase.component("c");
        auto const& dissolved_contam_component = liquid_phase.component("c");

        auto const porosity =
            medium.property(MaterialPropertyLib::PropertyType::porosity)
                .template value<double>(vars, pos, t, dt);

        auto const water_mol_mass =
            water_vapour_component
                .property(MaterialPropertyLib::PropertyType::molar_mass)
                .template value<double>(vars, pos, t, dt);
        auto const air_mol_mass =
            dry_air_component
                .property(MaterialPropertyLib::PropertyType::molar_mass)
                .template value<double>(vars, pos, t, dt);
        auto const contam_mol_mass =
            contam_vapour_component
                .property(MaterialPropertyLib::PropertyType::molar_mass)
                .template value<double>(vars, pos, t, dt);

        double const Sw =
            medium.property(MaterialPropertyLib::PropertyType::saturation)
                .template value<double>(vars, pos, t, dt);

        _saturation[ip] = Sw;
        vars[static_cast<int>(
            MaterialPropertyLib::Variable::liquid_saturation)] = Sw;

        double const dSw_dpc =
            medium.property(MaterialPropertyLib::PropertyType::saturation)
                .template dValue<double>(
                    vars, MaterialPropertyLib::Variable::capillary_pressure,
                    pos, t, dt);

        auto const density_water =
            liquid_phase.property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(vars, pos, t, dt);
        double const mol_density_water = density_water / water_mol_mass;
        double const mol_density_wet = mol_density_water;

        double const mol_density_nonwet =
            pg_int_pt / ideal_gas_constant_times_T_int_pt;
        double const mol_density_tot =
            Sw * mol_density_wet + (1 - Sw) * mol_density_nonwet;

        // specific latent heat of evaporation
        double const latent_heat_evaporation =
            water_vapour_component
                .property(
                    MaterialPropertyLib::PropertyType::specific_latent_heat)
                .template value<double>(vars, pos, t, dt);

        vars[static_cast<int>(
            MaterialPropertyLib::Variable::enthalpy_of_evaporation)] =
            latent_heat_evaporation;

        // saturated vapour pressure
        double const p_sat =
            water_vapour_component
                .property(MaterialPropertyLib::PropertyType::vapour_pressure)
                .template value<double>(vars, pos, t, dt);
        double const dp_sat_dT =
            water_vapour_component
                .property(MaterialPropertyLib::PropertyType::vapour_pressure)
                .template dValue<double>(
                    vars, MaterialPropertyLib::Variable::temperature, pos, t,
                    dt);

        // Kelvin-Laplace correction for menisci
        double const K = std::exp(-pc_int_pt / mol_density_water /
                                  ideal_gas_constant_times_T_int_pt);
        double const dK_dT = pc_int_pt / mol_density_water /
                             ideal_gas_constant_times_T_int_pt / T_int_pt * K;

        // vapour pressure inside pore space (water partial pressure in gas
        // phase)
        double const p_vapour_nonwet = p_sat * K;
        double const d_p_vapour_nonwet_dT = dp_sat_dT * K + p_sat * dK_dT;
        double const d_p_vapour_nonwet_dpc =
            p_vapour_nonwet *
            (-1 / mol_density_water / ideal_gas_constant_times_T_int_pt);

        // Henry constant of organic contaminant
        double const henry_contam =
            contam_vapour_component
                .property(MaterialPropertyLib::PropertyType::henry_constant)
                .template value<double>(vars, pos, t, dt);
        double d_henry_contam_dT =
            contam_vapour_component
                .property(MaterialPropertyLib::PropertyType::henry_constant)
                .template dValue<double>(
                    vars, MaterialPropertyLib::Variable::temperature, pos, t,
                    dt);

        // mass distribution coefficients of contam. and water
        double const k_c = pg_int_pt * henry_contam / mol_density_wet;
        double const k_w = pg_int_pt / p_vapour_nonwet;

        // intermediate parameter
        double const Ntot_c =
            Sw * mol_density_wet * k_c + (1 - Sw) * mol_density_nonwet;
        // phase-wise component molar fractions
        double const x_contam_nonwet = Xc_int_pt * mol_density_tot / Ntot_c;
        double const x_contam_wet = k_c * x_contam_nonwet;
        double const x_water_wet = 1 - x_contam_wet;
        double const x_water_nonwet = x_water_wet / k_w;
        double const x_air_nonwet = 1 - x_water_nonwet - x_contam_nonwet;

        _gas_molar_fraction_water[ip] = x_water_nonwet;
        _liquid_molar_fraction_contaminant[ip] = x_contam_wet;
        _gas_molar_fraction_contaminant[ip] = x_contam_nonwet;

        double const d_mol_density_nonwet_dpg =
            1 / ideal_gas_constant_times_T_int_pt;
        double const d_mol_density_nonwet_dT = -mol_density_nonwet / T_int_pt;
        double const d_mol_density_tot_dpc =
            (mol_density_wet - mol_density_nonwet) * dSw_dpc;
        double const d_mol_density_tot_dpg =
            (1 - Sw) * d_mol_density_nonwet_dpg;
        double const d_mol_density_tot_dT = (1 - Sw) * d_mol_density_nonwet_dT;

        double const d_kc_dpg = henry_contam / mol_density_wet;
        double const d_kc_dT = pg_int_pt * d_henry_contam_dT / mol_density_wet;

        double const d_x_contam_nonwet_dpc =
            Xc_int_pt * (d_mol_density_tot_dpc / Ntot_c -
                         (mol_density_wet * k_c - mol_density_nonwet) *
                             dSw_dpc * mol_density_tot / Ntot_c / Ntot_c);
        double const d_x_contam_nonwet_dpg =
            Xc_int_pt * (d_mol_density_tot_dpg / Ntot_c -
                         (Sw * mol_density_wet * d_kc_dpg +
                          (1 - Sw) * d_mol_density_nonwet_dpg) *
                             mol_density_tot / Ntot_c / Ntot_c);
        double const d_x_contam_nonwet_dXc = mol_density_tot / Ntot_c;
        double const d_x_contam_nonwet_dT =
            Xc_int_pt * (d_mol_density_tot_dT / Ntot_c -
                         (Sw * mol_density_wet * d_kc_dT +
                          (1 - Sw) * d_mol_density_nonwet_dT) *
                             mol_density_tot / Ntot_c / Ntot_c);

        double const d_x_contam_wet_dpc = k_c * d_x_contam_nonwet_dpc;
        double const d_x_contam_wet_dpg =
            k_c * d_x_contam_nonwet_dpg + d_kc_dpg * x_contam_nonwet;
        double const d_x_contam_wet_dXc = k_c * d_x_contam_nonwet_dXc;
        double const d_x_contam_wet_dT =
            k_c * d_x_contam_nonwet_dT + d_kc_dT * x_contam_nonwet;

        double const d_x_water_wet_dpc = -d_x_contam_wet_dpc;
        double const d_x_water_wet_dpg = -d_x_contam_wet_dpg;
        double const d_x_water_wet_dXc = -d_x_contam_wet_dXc;
        double const d_x_water_wet_dT = -d_x_contam_wet_dT;

        double const d_x_water_nonwet_dpc =
            (d_p_vapour_nonwet_dpc * x_water_wet +
             p_vapour_nonwet * d_x_water_wet_dpc) /
            pg_int_pt;
        double const d_x_water_nonwet_dpg =
            p_vapour_nonwet * (d_x_water_wet_dpg / pg_int_pt -
                               x_water_wet / pg_int_pt / pg_int_pt);
        double const d_x_water_nonwet_dXc = d_x_water_wet_dXc / k_w;
        double const d_x_water_nonwet_dT =
            (d_p_vapour_nonwet_dT * x_water_wet +
             p_vapour_nonwet * d_x_water_wet_dT) /
            pg_int_pt;

        double const d_x_air_nonwet_dpc =
            -d_x_water_nonwet_dpc - d_x_contam_nonwet_dpc;
        double const d_x_air_nonwet_dpg =
            -d_x_water_nonwet_dpg - d_x_contam_nonwet_dpg;
        double const d_x_air_nonwet_dXc =
            -d_x_water_nonwet_dXc - d_x_contam_nonwet_dXc;
        double const d_x_air_nonwet_dT =
            -d_x_water_nonwet_dT - d_x_contam_nonwet_dT;

        // mass densities
        double const density_contam_wet =
            mol_density_wet * contam_mol_mass * x_contam_wet;
        double const density_wet = density_water + density_contam_wet;

        double const density_air_nonwet =
            mol_density_nonwet * air_mol_mass * x_air_nonwet;
        double const density_water_nonwet =
            mol_density_nonwet * water_mol_mass * x_water_nonwet;
        double const density_contam_nonwet =
            mol_density_nonwet * contam_mol_mass * x_contam_nonwet;
        double const density_nonwet =
            density_air_nonwet + density_water_nonwet + density_contam_nonwet;

        auto const density_solid =
            solid_phase.property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(vars, pos, t, dt);

        // Derivative of nonwet phase density in terms of PVs
        double const d_density_nonwet_dpg =
            d_mol_density_nonwet_dpg *
                (air_mol_mass * x_air_nonwet + water_mol_mass * x_water_nonwet +
                 contam_mol_mass * x_contam_nonwet) +
            mol_density_nonwet * (air_mol_mass * d_x_air_nonwet_dpg +
                                  water_mol_mass * d_x_water_nonwet_dpg +
                                  contam_mol_mass * d_x_contam_nonwet_dpg);
        double const d_density_nonwet_dpc =
            mol_density_nonwet * (air_mol_mass * d_x_air_nonwet_dpc +
                                  water_mol_mass * d_x_water_nonwet_dpc +
                                  contam_mol_mass * d_x_contam_nonwet_dpc);
        double const d_density_nonwet_dXc =
            mol_density_nonwet * (air_mol_mass * d_x_air_nonwet_dXc +
                                  water_mol_mass * d_x_water_nonwet_dXc +
                                  contam_mol_mass * d_x_contam_nonwet_dXc);
        double const d_density_nonwet_dT =
            d_mol_density_nonwet_dT *
                (air_mol_mass * x_air_nonwet + water_mol_mass * x_water_nonwet +
                 contam_mol_mass * x_contam_nonwet) +
            mol_density_nonwet * (air_mol_mass * d_x_air_nonwet_dT +
                                  water_mol_mass * d_x_water_nonwet_dT +
                                  contam_mol_mass * d_x_contam_nonwet_dT);

        double const mol_mass_nonwet = x_water_nonwet * water_mol_mass +
                                       x_air_nonwet * air_mol_mass +
                                       x_contam_nonwet * contam_mol_mass;
        double const X_water_nonwet =
            x_water_nonwet * water_mol_mass / mol_mass_nonwet;
        double const X_air_nonwet =
            x_air_nonwet * air_mol_mass / mol_mass_nonwet;
        double const X_contam_nonwet = 1 - X_water_nonwet - X_air_nonwet;

        // heat capacity nonwetting phase
        double const heat_capacity_dry_air =
            dry_air_component
                .property(
                    MaterialPropertyLib::PropertyType::specific_heat_capacity)
                .template value<double>(vars, pos, t, dt);
        double const heat_capacity_water_vapour =
            water_vapour_component
                .property(
                    MaterialPropertyLib::PropertyType::specific_heat_capacity)
                .template value<double>(vars, pos, t, dt);
        double const heat_capacity_contam_vapour =
            contam_vapour_component
                .property(
                    MaterialPropertyLib::PropertyType::specific_heat_capacity)
                .template value<double>(vars, pos, t, dt);

        double const heat_capacity_water =
            liquid_phase
                .property(
                    MaterialPropertyLib::PropertyType::specific_heat_capacity)
                .template value<double>(vars, pos, t, dt);
        double const heat_capacity_solid =
            solid_phase
                .property(
                    MaterialPropertyLib::PropertyType::specific_heat_capacity)
                .template value<double>(vars, pos, t, dt);

        double const enthalpy_air_nonwet =
            heat_capacity_dry_air * (T_int_pt - CelsiusZeroInKelvin) +
            IdealGasConstant * (T_int_pt - CelsiusZeroInKelvin) / air_mol_mass;
        double const enthalpy_water_nonwet =
            heat_capacity_water_vapour * (T_int_pt - CelsiusZeroInKelvin) +
            latent_heat_evaporation;
        double const enthalpy_contam_nonwet =
            heat_capacity_contam_vapour * (T_int_pt - CelsiusZeroInKelvin) +
            IdealGasConstant * T_int_pt / contam_mol_mass;

        double const enthalpy_nonwet = enthalpy_air_nonwet * X_air_nonwet +
                                       enthalpy_water_nonwet * X_water_nonwet +
                                       enthalpy_contam_nonwet * X_contam_nonwet;
        double const enthalpy_wet =
            heat_capacity_water * (T_int_pt - CelsiusZeroInKelvin);

        double const internal_energy_nonwet =
            enthalpy_nonwet - pg_int_pt / density_nonwet;
        double const internal_energy_wet = enthalpy_wet;

        /// Derivatives
        double const d_enthalpy_air_nonwet_dT =
            heat_capacity_dry_air + IdealGasConstant / air_mol_mass;
        double const d_enthalpy_contaminant_nonwet_dT =
            heat_capacity_contam_vapour + IdealGasConstant / contam_mol_mass;

        double const d_enthalpy_nonwet_dT =
            heat_capacity_water * X_water_nonwet +
            d_enthalpy_air_nonwet_dT * X_air_nonwet +
            d_enthalpy_contaminant_nonwet_dT * X_contam_nonwet;

        // Assemble M matrix
        Map.noalias() += porosity * (1 - Sw) *
                         (mol_density_nonwet * d_x_air_nonwet_dpg +
                          x_air_nonwet * d_mol_density_nonwet_dpg) *
                         mass_operator;
        Mapc.noalias() +=
            porosity * mol_density_nonwet *
            ((1 - Sw) * d_x_air_nonwet_dpc - x_air_nonwet * dSw_dpc) *
            mass_operator;
        Max.noalias() += porosity * mol_density_nonwet * (1 - Sw) *
                         d_x_air_nonwet_dXc * mass_operator;
        Mat.noalias() += porosity * (1 - Sw) *
                         (mol_density_nonwet * d_x_air_nonwet_dT +
                          x_air_nonwet * d_mol_density_nonwet_dT) *
                         mass_operator;

        Mwp.noalias() +=
            porosity *
            (mol_density_wet * Sw * d_x_water_wet_dpg +
             (1 - Sw) * x_water_nonwet * d_mol_density_nonwet_dpg +
             mol_density_nonwet * (1 - Sw) * d_x_water_nonwet_dpg) *
            mass_operator;
        Mwpc.noalias() +=
            porosity *
            (mol_density_wet *
                 (x_water_wet * dSw_dpc + Sw * d_x_water_wet_dpc) +
             mol_density_nonwet *
                 ((1 - Sw) * d_x_water_nonwet_dpc - x_water_nonwet * dSw_dpc)) *
            mass_operator;
        Mwx.noalias() +=
            porosity *
            (mol_density_wet * Sw * d_x_water_wet_dXc +
             mol_density_nonwet * (1 - Sw) * d_x_water_nonwet_dXc) *
            mass_operator;
        Mwt.noalias() += porosity *
                         ((1 - Sw) / ideal_gas_constant_times_T_int_pt *
                          (d_p_vapour_nonwet_dT - p_vapour_nonwet / T_int_pt)) *
                         mass_operator;

        Mcp.noalias() +=
            porosity *
            (mol_density_wet * Sw * d_x_contam_wet_dpg +
             (1 - Sw) * x_contam_nonwet * d_mol_density_nonwet_dpg +
             mol_density_nonwet * (1 - Sw) * d_x_contam_nonwet_dpg) *
            mass_operator;
        Mcpc.noalias() +=
            porosity *
            (mol_density_wet *
                 (x_contam_wet * dSw_dpc + Sw * d_x_contam_wet_dpc) +
             mol_density_nonwet * ((1 - Sw) * d_x_contam_nonwet_dpc -
                                   x_contam_nonwet * dSw_dpc)) *
            mass_operator;
        Mcx.noalias() +=
            porosity *
            (mol_density_wet * Sw * d_x_contam_wet_dXc +
             mol_density_nonwet * (1 - Sw) * d_x_contam_nonwet_dXc) *
            mass_operator;
        Mct.noalias() +=
            porosity *
            (mol_density_wet * Sw * d_x_contam_wet_dT +
             (1 - Sw) * x_contam_nonwet * d_mol_density_nonwet_dT +
             mol_density_nonwet * (1 - Sw) * d_x_contam_nonwet_dT) *
            mass_operator;

        Mep.noalias() += porosity *
                         (d_density_nonwet_dpg * enthalpy_nonwet - 1) *
                         (1 - Sw) * mass_operator;
        Mepc.noalias() += porosity *
                              (density_wet * internal_energy_wet -
                               density_nonwet * internal_energy_nonwet) *
                              dSw_dpc * mass_operator +
                          porosity * d_density_nonwet_dpc * enthalpy_nonwet *
                              (1 - Sw) * mass_operator;
        Mex.noalias() += porosity * d_density_nonwet_dXc * enthalpy_nonwet *
                         (1 - Sw) * mass_operator;
        Met.noalias() +=
            ((1 - porosity) * density_solid * heat_capacity_solid +
             porosity * ((1 - Sw) * (d_density_nonwet_dT * enthalpy_nonwet +
                                     density_nonwet * d_enthalpy_nonwet_dT) +
                         Sw * density_wet * heat_capacity_water)) *
            mass_operator;

        double const diffusion_coeff_water_nonwet =
            water_vapour_component
                .property(MaterialPropertyLib::PropertyType::diffusion)
                .template value<double>(vars, pos, t, dt);
        double const diffusion_coeff_contam_nonwet =
            contam_vapour_component
                .property(MaterialPropertyLib::PropertyType::diffusion)
                .template value<double>(vars, pos, t, dt);
        double const diffusion_coeff_contam_wet =
            dissolved_contam_component
                .property(MaterialPropertyLib::PropertyType::diffusion)
                .template value<double>(vars, pos, t, dt);

        // nonwet
        double const k_rel_nonwet =
            medium
                .property(MaterialPropertyLib::PropertyType::
                              relative_permeability_nonwetting_phase)
                .template value<double>(vars, pos, t, dt);
        auto const mu_nonwet =
            gas_phase.property(MaterialPropertyLib::PropertyType::viscosity)
                .template value<double>(vars, pos, t, dt);
        double const lambda_nonwet = k_rel_nonwet / mu_nonwet;

        // wet
        double const k_rel_wet =
            medium
                .property(
                    MaterialPropertyLib::PropertyType::relative_permeability)
                .template value<double>(vars, pos, t, dt);
        auto const mu_wet =
            liquid_phase.property(MaterialPropertyLib::PropertyType::viscosity)
                .template value<double>(vars, pos, t, dt);
        double const lambda_wet = k_rel_wet / mu_wet;

        auto const permeability =
            MaterialPropertyLib::formEigenTensor<GlobalDim>(
                medium.property(MaterialPropertyLib::PropertyType::permeability)
                    .value(vars, pos, t, dt));

        auto const& b = _process_data.specific_body_force;

        GlobalDimVectorType const velocity_nonwet =
            _process_data.has_gravity
                ? GlobalDimVectorType(
                      -lambda_nonwet * permeability *
                      (dNdx * pg_nodal_values - density_nonwet * b))
                : GlobalDimVectorType(-lambda_nonwet * permeability * dNdx *
                                      pg_nodal_values);
        GlobalDimVectorType const velocity_wet =
            _process_data.has_gravity
                ? GlobalDimVectorType(
                      -lambda_wet * permeability *
                      (dNdx * (pg_nodal_values - pc_nodal_values) -
                       density_wet * b))
                : GlobalDimVectorType(-lambda_wet * permeability * dNdx *
                                      (pg_nodal_values - pc_nodal_values));

        laplace_operator.noalias() = dNdx.transpose() * permeability * dNdx * w;

        Ket.noalias() += w * N.transpose() *
                             (d_density_nonwet_dT * enthalpy_nonwet +
                              density_nonwet * d_enthalpy_nonwet_dT) *
                             velocity_nonwet.transpose() * dNdx +
                         w * N.transpose() * heat_capacity_water * density_wet *
                             velocity_wet.transpose() * dNdx;

        auto const solute_dispersivity_transverse =
            medium.template value<double>(
                MaterialPropertyLib::PropertyType::transversal_dispersivity);
        auto const solute_dispersivity_longitudinal =
            medium.template value<double>(
                MaterialPropertyLib::PropertyType::longitudinal_dispersivity);

        double const velocity_wet_magnitude = velocity_wet.norm();
        GlobalDimMatrixType const hydrodynamic_dispersion =
            velocity_wet_magnitude != 0.0
                ? GlobalDimMatrixType(
                      (porosity * Sw * diffusion_coeff_contam_wet +
                       solute_dispersivity_transverse *
                           velocity_wet_magnitude) *
                          I +
                      (solute_dispersivity_longitudinal -
                       solute_dispersivity_transverse) /
                          velocity_wet_magnitude * velocity_wet *
                          velocity_wet.transpose())
                : GlobalDimMatrixType(
                      (porosity * Sw * diffusion_coeff_contam_wet +
                       solute_dispersivity_transverse *
                           velocity_wet_magnitude) *
                      I);

        auto const dispersion_operator =
            dNdx.transpose() * hydrodynamic_dispersion * dNdx * w;

        // The sum of all diffusive fluxes in either phase must equal to zero
        Kap.noalias() +=
            (mol_density_nonwet * x_air_nonwet * lambda_nonwet) *
                laplace_operator -
            (1 - Sw) * porosity * mol_density_nonwet *
                (diffusion_coeff_water_nonwet * d_x_water_nonwet_dpg +
                 diffusion_coeff_contam_nonwet * d_x_contam_nonwet_dpg) *
                diffusion_operator;
        Kapc.noalias() +=
            -(1 - Sw) * porosity * mol_density_nonwet *
            (diffusion_coeff_water_nonwet * d_x_water_nonwet_dpc +
             diffusion_coeff_contam_nonwet * d_x_contam_nonwet_dpc) *
            diffusion_operator;
        Kax.noalias() +=
            -(1 - Sw) * porosity * mol_density_nonwet *
            (diffusion_coeff_water_nonwet * d_x_water_nonwet_dXc +
             diffusion_coeff_contam_nonwet * d_x_contam_nonwet_dXc) *
            diffusion_operator;
        Kat.noalias() +=
            -(1 - Sw) * porosity * mol_density_nonwet *
            (diffusion_coeff_water_nonwet * d_x_water_nonwet_dT +
             diffusion_coeff_contam_nonwet * d_x_contam_nonwet_dT) *
            diffusion_operator;

        Kwp.noalias() +=
            (mol_density_wet * x_water_wet * lambda_wet +
             mol_density_nonwet * x_water_nonwet * lambda_nonwet) *
                laplace_operator +
            porosity *
                (Sw * mol_density_wet * diffusion_coeff_contam_wet *
                     d_x_water_wet_dpg +
                 (1 - Sw) * diffusion_coeff_water_nonwet * mol_density_nonwet *
                     d_x_water_nonwet_dpg) *
                diffusion_operator;
        Kwpc.noalias() +=
            -mol_density_wet * x_water_wet * lambda_wet * laplace_operator +
            porosity *
                (Sw * mol_density_wet * diffusion_coeff_contam_wet *
                     d_x_water_wet_dpc +
                 (1 - Sw) * mol_density_nonwet * diffusion_coeff_water_nonwet *
                     d_x_water_nonwet_dpc) *
                diffusion_operator;
        Kwx.noalias() +=
            porosity *
            (Sw * mol_density_wet * diffusion_coeff_contam_wet *
                 d_x_water_wet_dXc +
             (1 - Sw) * mol_density_nonwet * diffusion_coeff_water_nonwet *
                 d_x_water_nonwet_dXc) *
            diffusion_operator;
        Kwt.noalias() +=
            porosity *
            (Sw * mol_density_wet * diffusion_coeff_contam_wet *
                 d_x_water_wet_dT +
             (1 - Sw) * mol_density_nonwet * diffusion_coeff_water_nonwet *
                 d_x_water_nonwet_dT) *
            diffusion_operator;

        Kcp.noalias() +=
            (mol_density_wet * x_contam_wet * lambda_wet +
             mol_density_nonwet * x_contam_nonwet * lambda_nonwet) *
                laplace_operator +
            mol_density_wet * dispersion_operator * d_x_contam_wet_dpg +
            porosity * (1 - Sw) * mol_density_nonwet *
                diffusion_coeff_contam_nonwet * d_x_contam_nonwet_dpg *
                diffusion_operator;
        Kcpc.noalias() +=
            -mol_density_wet * x_contam_wet * lambda_wet * laplace_operator +
            mol_density_wet * dispersion_operator * d_x_contam_wet_dpc +
            porosity * (1 - Sw) * mol_density_nonwet *
                diffusion_coeff_contam_nonwet * d_x_contam_nonwet_dpc *
                diffusion_operator;
        Kcx.noalias() +=
            mol_density_wet * dispersion_operator * d_x_contam_wet_dXc +
            porosity * (1 - Sw) * mol_density_nonwet *
                diffusion_coeff_contam_nonwet * d_x_contam_nonwet_dXc *
                diffusion_operator;
        Kct.noalias() +=
            mol_density_wet * dispersion_operator * d_x_contam_wet_dT +
            porosity * (1 - Sw) * mol_density_nonwet *
                diffusion_coeff_contam_nonwet * d_x_contam_nonwet_dT *
                diffusion_operator;

        Kep.noalias() += (lambda_nonwet * density_nonwet * enthalpy_nonwet +
                          lambda_wet * density_wet * enthalpy_wet) *
                             laplace_operator +
                         (1 - Sw) * porosity * diffusion_coeff_water_nonwet *
                             mol_density_nonwet *
                             (air_mol_mass * enthalpy_air_nonwet -
                              water_mol_mass * enthalpy_water_nonwet) *
                             d_x_air_nonwet_dpg * diffusion_operator;
        Kepc.noalias() +=
            -lambda_wet * enthalpy_wet * density_wet * laplace_operator +
            (1 - Sw) * porosity * diffusion_coeff_water_nonwet *
                mol_density_nonwet *
                (air_mol_mass * enthalpy_air_nonwet -
                 water_mol_mass * enthalpy_water_nonwet) *
                d_x_air_nonwet_dpc * diffusion_operator;

        if (medium.hasProperty(
                MaterialPropertyLib::PropertyType::thermal_conductivity))
        {
            auto const lambda =
                medium
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .value(vars, pos, t, dt);

            GlobalDimMatrixType const heat_conductivity_unsaturated =
                MaterialPropertyLib::formEigenTensor<GlobalDim>(lambda);

            Ket.noalias() +=
                dNdx.transpose() * heat_conductivity_unsaturated * dNdx * w +
                (1 - Sw) * porosity * diffusion_coeff_water_nonwet *
                    mol_density_nonwet *
                    (air_mol_mass * enthalpy_air_nonwet -
                     water_mol_mass * enthalpy_water_nonwet) *
                    d_x_air_nonwet_dT * diffusion_operator;
        }
        else
        {
            auto const thermal_conductivity_solid =
                solid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .value(vars, pos, t, dt);

            auto const thermal_conductivity_fluid =
                liquid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .template value<double>(vars, pos, t, dt) *
                Sw;

            GlobalDimMatrixType const heat_conductivity_unsaturated =
                MaterialPropertyLib::formEffectiveThermalConductivity<
                    GlobalDim>(thermal_conductivity_solid,
                               thermal_conductivity_fluid, porosity);

            Ket.noalias() +=
                dNdx.transpose() * heat_conductivity_unsaturated * dNdx * w +
                (1 - Sw) * porosity * diffusion_coeff_water_nonwet *
                    mol_density_nonwet *
                    (air_mol_mass * enthalpy_air_nonwet -
                     water_mol_mass * enthalpy_water_nonwet) *
                    d_x_air_nonwet_dT * diffusion_operator;
        }

        if (_process_data.has_gravity)
        {
            NodalVectorType gravity_operator =
                dNdx.transpose() * permeability * b * w;
            Ba.noalias() += (mol_density_nonwet * x_air_nonwet * lambda_nonwet *
                             density_nonwet) *
                            gravity_operator;
            Bw.noalias() +=
                (mol_density_wet * x_water_wet * lambda_wet * density_wet +
                 mol_density_nonwet * x_water_nonwet * lambda_nonwet *
                     density_nonwet) *
                gravity_operator;
            Bc.noalias() +=
                (mol_density_wet * x_contam_wet * lambda_wet * density_wet +
                 mol_density_nonwet * x_contam_nonwet * lambda_nonwet *
                     density_nonwet) *
                gravity_operator;
            Be.noalias() +=
                (lambda_nonwet * density_nonwet * density_nonwet *
                     enthalpy_nonwet +
                 lambda_wet * density_wet * density_wet * enthalpy_wet) *
                gravity_operator;
        }  // end of has gravity
    }
    if (_process_data.has_mass_lumping)
    {
        Map = Map.colwise().sum().eval().asDiagonal();
        Mapc = Mapc.colwise().sum().eval().asDiagonal();
        Max = Max.colwise().sum().eval().asDiagonal();
        Mat = Mat.colwise().sum().eval().asDiagonal();
        Mwp = Mwp.colwise().sum().eval().asDiagonal();
        Mwpc = Mwpc.colwise().sum().eval().asDiagonal();
        Mwx = Mwx.colwise().sum().eval().asDiagonal();
        Mwt = Mwt.colwise().sum().eval().asDiagonal();
        Mcp = Mcp.colwise().sum().eval().asDiagonal();
        Mcpc = Mcpc.colwise().sum().eval().asDiagonal();
        Mcx = Mcx.colwise().sum().eval().asDiagonal();
        Mct = Mct.colwise().sum().eval().asDiagonal();
        Mep = Mep.colwise().sum().eval().asDiagonal();
        Mepc = Mepc.colwise().sum().eval().asDiagonal();
        Mex = Mex.colwise().sum().eval().asDiagonal();
        Met = Met.colwise().sum().eval().asDiagonal();
    }  // end of mass-lumping
}

}  // namespace ThermalTwoPhaseFlowWithPP
}  // namespace ProcessLib
