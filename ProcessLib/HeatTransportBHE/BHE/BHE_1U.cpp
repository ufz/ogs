/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BHE_1U.h"

using namespace ProcessLib::HeatTransportBHE::BHE;

void BHE_1U::initialize()

{
    double tmp_u = calcPipeFlowVelocity(Q_r, 2.0 * pipe_param.r_inner);
    _u(0) = tmp_u;
    _u(1) = tmp_u;

    // calculate Renolds number
    Re = calcRenoldsNumber(std::abs(_u(0)),
                           2.0 * pipe_param.r_inner,
                           refrigerant_param.mu_r,
                           refrigerant_param.rho_r);
    // calculate Prandtl number
    Pr = calcPrandtlNumber(refrigerant_param.mu_r,
                           refrigerant_param.heat_cap_r,
                           refrigerant_param.lambda_r);

    // calculate Nusselt number
    double tmp_Nu =
        calcNusseltNumber(2.0 * pipe_param.r_inner, borehole_geometry.length);
    _Nu(0) = tmp_Nu;
    _Nu(1) = tmp_Nu;

    calcThermalResistances();
    calcHeatTransferCoefficients();
}

/**
 * calculate thermal resistance
 */
void BHE_1U::calcThermalResistances()
{
    // thermal resistance due to the grout transition
    double chi;
    double d0;  // the average outer diameter of the pipes
    // double s; // diagonal distances of pipes
    double R_adv, R_con;
    double const& D = borehole_geometry.diameter;
    // double const& L = borehole_geometry.L;
    double const& lambda_r = refrigerant_param.lambda_r;
    double const& lambda_g = grout_param.lambda_g;
    double const& lambda_p = pipe_param.lambda_p;

    // thermal resistance due to thermal conductivity of the pip wall material
    // Eq. 36 in Diersch_2011_CG
    _R_adv_i1 = 1.0 / (_Nu(0) * lambda_r * PI);
    _R_adv_o1 = 1.0 / (_Nu(1) * lambda_r * PI);

    d0 = 2.0 * pipe_param.r_outer;
    // s = omega * std::sqrt(2);
    // Eq. 49
    _R_con_a_i1 = _R_con_a_o1 =
        std::log(pipe_param.r_outer / pipe_param.r_inner) /
        (2.0 * PI * lambda_p);
    // Eq. 51
    chi = std::log(std::sqrt(D * D + 2 * d0 * d0) / 2 / d0) /
          std::log(D / std::sqrt(2) / d0);
    if (extern_Ra_Rb.use_extern_Ra_Rb)
    {
        R_adv = 0.5 * (_R_adv_i1 + _R_adv_o1);
        R_con = 0.5 * (_R_con_a_i1 + _R_con_a_o1);
        _R_g = 2 * extern_Ra_Rb.ext_Rb - R_adv - R_con;
    }
    else
    {
        // Eq. 52
        _R_g = acosh((D * D + d0 * d0 - omega * omega) / (2 * D * d0)) /
               (2 * PI * lambda_g) * (1.601 - 0.888 * omega / D);
    }
    _R_con_b = chi * _R_g;
    // Eq. 29 and 30
    if (extern_def_thermal_resistances.if_use_defined_therm_resis)
    {
        _R_fig = extern_def_thermal_resistances.ext_Rfig;
        _R_fog = extern_def_thermal_resistances.ext_Rfog;
    }
    else
    {
        _R_fig = _R_adv_i1 + _R_con_a_i1 + _R_con_b;
        _R_fog = _R_adv_o1 + _R_con_a_o1 + _R_con_b;
    }

    // thermal resistance due to grout-soil exchange
    if (extern_def_thermal_resistances.if_use_defined_therm_resis)
        _R_gs = extern_def_thermal_resistances.ext_Rgs;
    else
        _R_gs = (1 - chi) * _R_g;

    // thermal resistance due to inter-grout exchange
    double R_ar;
    if (extern_Ra_Rb.use_extern_Ra_Rb)
    {
        R_ar = extern_Ra_Rb.ext_Ra - 2 * (R_adv + R_con);
    }
    else
    {
        R_ar = acosh((2.0 * omega * omega - d0 * d0) / d0 / d0) /
               (2.0 * PI * lambda_g);
    }

    if (extern_def_thermal_resistances.if_use_defined_therm_resis)
        _R_gg = extern_def_thermal_resistances.ext_Rgg1;
    else
        _R_gg = 2.0 * _R_gs * (R_ar - 2.0 * chi * _R_g) /
                (2.0 * _R_gs - R_ar + 2.0 * chi * _R_g);

    if (!std::isfinite(_R_gg))
    {
        OGS_FATAL(
            "Error!!! Grout Thermal Resistance is an infinite number! The "
            "simulation will be stopped!");
    }

    // check if constraints regarding negative thermal resistances are violated
    // apply correction procedure
    // Section (1.5.5) in FEFLOW White Papers Vol V.
    double constraint = 1.0 / ((1.0 / _R_gg) + (1.0 / (2.0 * _R_gs)));
    int count = 0;
    while (constraint < 0.0)
    {
        if (extern_def_thermal_resistances.if_use_defined_therm_resis ||
            extern_Ra_Rb.use_extern_Ra_Rb)
        {
            OGS_FATAL(
                "Error!!! Constraints on thermal resistances are violated! "
                "Correction procedure can't be applied due to user defined "
                "thermal resistances! The simulation will be stopped!");
        }
        if (count == 0)
        {
            chi *= 0.66;
            _R_gs = (1 - chi) * _R_g;
            _R_gg = 2.0 * _R_gs * (R_ar - 2.0 * chi * _R_g) /
                    (2.0 * _R_gs - R_ar + 2.0 * chi * _R_g);
        }
        if (count == 1)
        {
            chi *= 0.5;
            _R_gs = (1 - chi) * _R_g;
            _R_gg = 2.0 * _R_gs * (R_ar - 2.0 * chi * _R_g) /
                    (2.0 * _R_gs - R_ar + 2.0 * chi * _R_g);
        }
        if (count == 2)
        {
            chi = 0.0;
            _R_gs = (1 - chi) * _R_g;
            _R_gg = 2.0 * _R_gs * (R_ar - 2.0 * chi * _R_g) /
                    (2.0 * _R_gs - R_ar + 2.0 * chi * _R_g);
            break;
        }
        DBUG(
            "Warning! Correction procedure was applied due to negative thermal "
            "resistance! Correction step #%d.\n",
            count);
        constraint = 1.0 / ((1.0 / _R_gg) + (1.0 / (2.0 * _R_gs)));
        count++;
    }

    // kleep the following lines------------------------------------------------
    // when debugging the code, printing the R and phi values are needed--------
    // std::cout << "Rfig =" << _R_fig << " Rfog =" << _R_fog << " Rgg =" <<
    // _R_gg << " Rgs =" << _R_gs << "\n"; double phi_fig = 1.0 / (_R_fig *
    // S_i); double phi_fog = 1.0 / (_R_fog * S_o); double phi_gg = 1.0 / (_R_gg
    // * S_g1); double phi_gs = 1.0 / (_R_gs * S_gs); std::cout << "phi_fig ="
    // << phi_fig << " phi_fog =" << phi_fog << " phi_gg =" << phi_gg << "
    // phi_gs =" << phi_gs << "\n";
    // -------------------------------------------------------------------------
}

/**
 * calculate heat transfer coefficient
 */
void BHE_1U::calcHeatTransferCoefficients()
{
    _PHI_fig = 1.0 / _R_fig;
    _PHI_fog = 1.0 / _R_fog;
    _PHI_gg = 1.0 / _R_gg;
    _PHI_gs = 1.0 / _R_gs;
}

double BHE_1U::getMassCoeff(std::size_t idx_unknown) const
{
    double const& rho_r = refrigerant_param.rho_r;
    double const& heat_cap_r = refrigerant_param.heat_cap_r;
    double const& rho_g = grout_param.rho_g;
    double const& porosity_g = grout_param.porosity_g;
    double const& heat_cap_g = grout_param.heat_cap_g;

    double mass_coeff = 0.0;

    switch (idx_unknown)
    {
        case 0:  // i1
            mass_coeff = rho_r * heat_cap_r * CSA_i;
            break;
        case 1:  // o1
            mass_coeff = rho_r * heat_cap_r * CSA_o;
            break;
        case 2:  // g1
            mass_coeff = (1.0 - porosity_g) * rho_g * heat_cap_g * CSA_g1;
            break;
        case 3:  // g2
            mass_coeff = (1.0 - porosity_g) * rho_g * heat_cap_g * CSA_g2;
            break;
        default:
            break;
    }

    return mass_coeff;
}

void BHE_1U::getLaplaceMatrix(std::size_t idx_unknown,
                              Eigen::MatrixXd& mat_laplace) const
{
    double const& lambda_r = refrigerant_param.lambda_r;
    double const& rho_r = refrigerant_param.rho_r;
    double const& heat_cap_r = refrigerant_param.heat_cap_r;
    double const& alpha_L = refrigerant_param.alpha_L;
    double const& porosity_g = grout_param.porosity_g;
    double const& lambda_g = grout_param.lambda_g;

    // Here we calculates the laplace coefficients in the governing
    // equations of BHE. These governing equations can be found in
    // 1) Diersch (2013) FEFLOW book on page 952, M.120-122, or
    // 2) Diersch (2011) Comp & Geosci 37:1122-1135, Eq. 19-22.
    double laplace_coeff(0.0);
    mat_laplace.setZero();

    switch (idx_unknown)
    {
        case 0:
            // pipe i1, Eq. 19
            laplace_coeff =
                (lambda_r + rho_r * heat_cap_r * alpha_L * _u.norm()) * CSA_i;
            break;
        case 1:
            // pipe o1, Eq. 20
            laplace_coeff =
                (lambda_r + rho_r * heat_cap_r * alpha_L * _u.norm()) * CSA_o;
            break;
        case 2:
            // pipe g1, Eq. 21
            laplace_coeff = (1.0 - porosity_g) * lambda_g * CSA_g1;
            break;
        case 3:
            // pipe g2, Eq. 22
            laplace_coeff = (1.0 - porosity_g) * lambda_g * CSA_g2;
            break;
        default:
            OGS_FATAL(
                "Error !!! The index passed to get_laplace_coeff for BHE is "
                "not correct. ");
            break;
    }

    mat_laplace(0, 0) = laplace_coeff;
    mat_laplace(1, 1) = laplace_coeff;
    mat_laplace(2, 2) = laplace_coeff;
}

void BHE_1U::getAdvectionVector(std::size_t idx_unknown,
                                Eigen::VectorXd& vec_advection) const
{
    double const& rho_r = refrigerant_param.rho_r;
    double const& heat_cap_r = refrigerant_param.heat_cap_r;
    double advection_coeff(0);
    vec_advection.setZero();

    switch (idx_unknown)
    {
        case 0:
            // pipe i1, Eq. 19
            advection_coeff = rho_r * heat_cap_r * _u(0) * CSA_i;
            // z direction
            vec_advection(2) = -1.0 * advection_coeff;
            break;
        case 1:
            // pipe o1, Eq. 20
            advection_coeff = rho_r * heat_cap_r * _u(0) * CSA_o;
            // z direction
            vec_advection(2) = 1.0 * advection_coeff;
            break;
        case 2:
            // grout g1, Eq. 21
            advection_coeff = 0.0;
            break;
        case 3:
            // grout g2, Eq. 22
            advection_coeff = 0.0;
            break;
        default:
            OGS_FATAL(
                "Error !!! The index passed to get_advection_coeff for BHE is "
                "not correct. ");
            break;
    }
}

void ProcessLib::HeatTransportBHE::BHE::BHE_1U::setRMatrices(
    const int idx_bhe_unknowns, const int NumNodes,
    Eigen::MatrixXd& matBHE_loc_R,
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>&
        R_matrix,
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>&
        R_pi_s_matrix,
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>&
        R_s_matrix) const
{
    switch (idx_bhe_unknowns)
    {
        case 0:  // PHI_fig
            R_matrix.block(0, 2 * NumNodes, NumNodes, NumNodes) +=
                -1.0 * matBHE_loc_R;
            R_matrix.block(2 * NumNodes, 0, NumNodes, NumNodes) +=
                -1.0 * matBHE_loc_R;

            R_matrix.block(0, 0, NumNodes, NumNodes) +=
                1.0 * matBHE_loc_R;  // K_i1
            R_matrix.block(2 * NumNodes, 2 * NumNodes, NumNodes, NumNodes) +=
                1.0 * matBHE_loc_R;  // K_ig
            break;
        case 1:  // PHI_fog
            R_matrix.block(NumNodes, 3 * NumNodes, NumNodes, NumNodes) +=
                -1.0 * matBHE_loc_R;
            R_matrix.block(3 * NumNodes, NumNodes, NumNodes, NumNodes) +=
                -1.0 * matBHE_loc_R;

            R_matrix.block(NumNodes, NumNodes, NumNodes, NumNodes) +=
                1.0 * matBHE_loc_R;  // K_o1
            R_matrix.block(3 * NumNodes, 3 * NumNodes, NumNodes, NumNodes) +=
                1.0 * matBHE_loc_R;  // K_og
            break;
        case 2:  // PHI_gg
            R_matrix.block(2 * NumNodes, 3 * NumNodes, NumNodes, NumNodes) +=
                -1.0 * matBHE_loc_R;
            R_matrix.block(3 * NumNodes, 2 * NumNodes, NumNodes, NumNodes) +=
                -1.0 * matBHE_loc_R;

            R_matrix.block(2 * NumNodes, 2 * NumNodes, NumNodes, NumNodes) +=
                1.0 * matBHE_loc_R;  // K_ig  // notice we only have
                                     // 1 PHI_gg term here.
            R_matrix.block(3 * NumNodes, 3 * NumNodes, NumNodes, NumNodes) +=
                1.0 * matBHE_loc_R;  // K_og  // see Diersch 2013 FEFLOW
                                     // book page 954 Table M.2
            break;
        case 3:  // PHI_gs
            R_s_matrix += 1.0 * matBHE_loc_R;

            R_pi_s_matrix.block(2 * NumNodes, 0, NumNodes, NumNodes) +=
                -1.0 * matBHE_loc_R;
            R_pi_s_matrix.block(3 * NumNodes, 0, NumNodes, NumNodes) +=
                -1.0 * matBHE_loc_R;
            R_matrix.block(2 * NumNodes, 2 * NumNodes, NumNodes, NumNodes) +=
                1.0 * matBHE_loc_R;  // K_ig
            R_matrix.block(3 * NumNodes, 3 * NumNodes, NumNodes, NumNodes) +=
                1.0 * matBHE_loc_R;  // K_og
            break;
    }
}

double BHE_1U::getBoundaryHeatExchangeCoeff(std::size_t idx_unknown) const
{
    // Here we calculates the boundary heat exchange coefficients
    // in the governing equations of BHE.
    // These governing equations can be found in
    // 1) Diersch (2013) FEFLOW book on page 958, M.3, or
    // 2) Diersch (2011) Comp & Geosci 37:1122-1135, Eq. 90-97.

    double exchange_coeff(0);

    switch (idx_unknown)
    {
        case 0:
            // PHI_fig
            exchange_coeff = _PHI_fig;
            break;
        case 1:
            // PHI_fog
            exchange_coeff = _PHI_fog;
            break;
        case 2:
            // PHI_gg
            exchange_coeff = _PHI_gg;
            break;
        case 3:
            // PHI_gs
            exchange_coeff = _PHI_gs;
            break;
        default:
            OGS_FATAL(
                "Error !!! The index passed to "
                "getBoundaryHeatExchangeCoeff for BHE is not correct. ");
            break;
    }
    return exchange_coeff;
}

double BHE_1U::getTinByTout(double T_out, double current_time = -1.0)
{
    double T_in(0.0);
    double power_tmp(0.0);
    double building_power_tmp(0.0);
    // double power_elect_tmp(0.0);
    double Q_r_tmp(0.0);
    double COP_tmp(0.0);
    double fac_dT = 1.0;
    double const& rho_r = refrigerant_param.rho_r;
    double const& heat_cap_r = refrigerant_param.heat_cap_r;

    switch (this->boundary_type)
    {
        case BHE_BOUNDARY_TYPE::FIXED_INFLOW_TEMP_CURVE_BOUNDARY:
        {
            T_in = inflow_temperature_curve->getValue(current_time);
        }
        break;
        case BHE_BOUNDARY_TYPE::POWER_IN_WATT_BOUNDARY:
            if (use_flowrate_curve)
            {
                Q_r_tmp = flowrate_curve->getValue(current_time);

                updateFlowRate(Q_r_tmp);
            }
            else
                Q_r_tmp = Q_r;
            T_in = power_in_watt_val / Q_r_tmp / heat_cap_r / rho_r + T_out;
            break;
        case BHE_BOUNDARY_TYPE::FIXED_TEMP_DIFF_BOUNDARY:
            if (use_flowrate_curve)
            {
                Q_r_tmp = flowrate_curve->getValue(current_time);

                updateFlowRate(Q_r_tmp);
            }
            T_in = T_out + delta_T_val;
            break;
        case BHE_BOUNDARY_TYPE::POWER_IN_WATT_CURVE_FIXED_DT_BOUNDARY:
            // get the power value in the curve
            // power_tmp = GetCurveValue(power_in_watt_curve_idx, 0,
            // current_time, &flag_valid);
            power_tmp = power_in_watt_curve->getValue(current_time);

            if (power_tmp < 0)
                fac_dT = -1.0;
            else
                fac_dT = 1.0;
            // if power value exceeds threshold, calculate new values
            if (fabs(power_tmp) > threshold)
            {
                // calculate the corresponding flow rate needed
                // using the defined delta_T value
                Q_r_tmp =
                    power_tmp / (fac_dT * delta_T_val) / heat_cap_r / rho_r;
                // update all values dependent on the flow rate
                updateFlowRate(Q_r_tmp);
                // calculate the new T_in
                T_in = T_out + (fac_dT * delta_T_val);
                // print out updated flow rate
                // std::cout << "Qr: " << Q_r_tmp << std::endl;
            }
            else
            {
                Q_r_tmp = 1.0e-12;  // this has to be a small value to avoid
                                    // division by zero
                // update all values dependent on the flow rate
                updateFlowRate(Q_r_tmp);
                // calculate the new T_in
                T_in = T_out;
                // print out updated flow rate
                // std::cout << "Qr: " << Q_r_tmp << std::endl;
            }
            break;
        case BHE_BOUNDARY_TYPE::BUILDING_POWER_IN_WATT_CURVE_FIXED_DT_BOUNDARY:
            // get the building power value in the curve
            // building_power_tmp = GetCurveValue(power_in_watt_curve_idx, 0,
            // current_time, &flag_valid);
            building_power_tmp = power_in_watt_curve->getValue(current_time);

            if (building_power_tmp <= 0.0)
            {
                // get COP value based on T_out in the curve
                COP_tmp = heating_cop_curve->getValue(T_out);

                // now calculate how much power needed from BHE
                power_tmp = building_power_tmp * (COP_tmp - 1.0) / COP_tmp;
                // also how much power from electricity
                // power_elect_tmp = building_power_tmp - power_tmp;
                // print the amount of power needed
                // std::cout << "COP: " << COP_tmp << ", Q_bhe: " << power_tmp
                // << ", Q_elect: " << power_elect_tmp << std::endl;
                fac_dT = -1.0;
            }
            else
            {
                // get COP value based on T_out in the curve
                COP_tmp = cooling_cop_curve->getValue(T_out);

                // now calculate how much power needed from BHE
                power_tmp = building_power_tmp * (COP_tmp + 1.0) / COP_tmp;
                // also how much power from electricity
                // power_elect_tmp = -building_power_tmp + power_tmp;
                // print the amount of power needed
                // std::cout << "COP: " << COP_tmp << ", Q_bhe: " << power_tmp
                // << ", Q_elect: " << power_elect_tmp << std::endl;
                fac_dT = 1.0;
            }
            // if power value exceeds threshold, calculate new values
            if (fabs(power_tmp) > threshold)
            {
                // calculate the corresponding flow rate needed
                // using the defined delta_T value
                Q_r_tmp =
                    power_tmp / (fac_dT * delta_T_val) / heat_cap_r / rho_r;
                // update all values dependent on the flow rate
                updateFlowRate(Q_r_tmp);
                // calculate the new T_in
                T_in = T_out + (fac_dT * delta_T_val);
                // print out updated flow rate
                // std::cout << "Qr: " << Q_r_tmp << std::endl;
            }
            else
            {
                Q_r_tmp = 1.0e-12;  // this has to be a small value to avoid
                                    // division by zero
                // update all values dependent on the flow rate
                updateFlowRate(Q_r_tmp);
                // calculate the new T_in
                T_in = T_out;
                // print out updated flow rate
                // std::cout << "Qr: " << Q_r_tmp << std::endl;
            }
            break;
        case BHE_BOUNDARY_TYPE::
            BUILDING_POWER_IN_WATT_CURVE_FIXED_FLOW_RATE_BOUNDARY:
            // get the building power value in the curve
            // building_power_tmp = GetCurveValue(power_in_watt_curve_idx, 0,
            // current_time, &flag_valid);
            building_power_tmp = power_in_watt_curve->getValue(current_time);

            if (building_power_tmp <= 0)
            {
                // get COP value based on T_out in the curve
                COP_tmp = heating_cop_curve->getValue(T_out);
                // now calculate how much power needed from BHE
                power_tmp = building_power_tmp * (COP_tmp - 1.0) / COP_tmp;
                // also how much power from electricity
                // power_elect_tmp = building_power_tmp - power_tmp;
                // print the amount of power needed
                // std::cout << "COP: " << COP_tmp << ", Q_bhe: " << power_tmp
                // << ", Q_elect: " << power_elect_tmp << std::endl;
            }
            else
            {
                // get COP value based on T_out in the curve
                COP_tmp = cooling_cop_curve->getValue(T_out);
                // now calculate how much power needed from BHE
                power_tmp = building_power_tmp * (COP_tmp + 1.0) / COP_tmp;
                // also how much power from electricity
                // power_elect_tmp = -building_power_tmp + power_tmp;
                // print the amount of power needed
                // std::cout << "COP: " << COP_tmp << ", Q_bhe: " << power_tmp
                // << ", Q_elect: " << power_elect_tmp << std::endl;
            }
            // Assign Qr whether from curve or fixed value
            if (use_flowrate_curve)
            {
                // Q_r_tmp = GetCurveValue(flowrate_curve_idx, 0, current_time,
                // &flag_valid);
                Q_r_tmp = flowrate_curve->getValue(current_time);
                updateFlowRate(Q_r_tmp);
            }
            else
                Q_r_tmp = Q_r;
            if (fabs(power_tmp) < threshold)
            {
                Q_r_tmp = 1.0e-12;  // this has to be a small value to avoid
                                    // division by zero update all values
                                    // dependent on the flow rate
                updateFlowRate(Q_r_tmp);
                // calculate the new T_in
                T_in = T_out;
                // print out updated flow rate
                // std::cout << "Qr: " << Q_r_tmp << std::endl;
            }
            else
            {
                Q_r_tmp = Q_r;
                updateFlowRate(Q_r_tmp);
                // calculate the dT value based on fixed flow rate
                delta_T_val = power_tmp / Q_r_tmp / heat_cap_r / rho_r;
                // calcuate the new T_in
                T_in = T_out + delta_T_val;
            }
            break;
        case BHE_BOUNDARY_TYPE::POWER_IN_WATT_CURVE_FIXED_FLOW_RATE_BOUNDARY:
            // get the power value in the curve
            // power_tmp = GetCurveValue(power_in_watt_curve_idx, 0,
            // current_time, &flag_valid);
            power_tmp = power_in_watt_curve->getValue(current_time);

            // Assign Qr whether from curve or fixed value
            if (use_flowrate_curve)
            {
                // Q_r_tmp = GetCurveValue(flowrate_curve_idx, 0, current_time,
                // &flag_valid);
                Q_r_tmp = flowrate_curve->getValue(current_time);
                updateFlowRate(Q_r_tmp);
            }
            else
                Q_r_tmp = Q_r;
            // calculate the dT value based on fixed flow rate
            if (fabs(power_tmp) < threshold)
            {
                Q_r_tmp = 1.0e-12;  // this has to be a small value to avoid
                                    // division by zero update all values
                                    // dependent on the flow rate
                updateFlowRate(Q_r_tmp);
                // calculate the new T_in
                T_in = T_out;
                // print out updated flow rate
                // std::cout << "Qr: " << Q_r_tmp << std::endl;
            }
            else
            {
                Q_r_tmp = Q_r;
                updateFlowRate(Q_r_tmp);
                // calculate the dT value based on fixed flow rate
                delta_T_val = power_tmp / Q_r_tmp / heat_cap_r / rho_r;
                // calcuate the new T_in
                T_in = T_out + delta_T_val;
            }
            break;
        default:
            T_in = T_out;
            break;
    }

    return T_in;
}
