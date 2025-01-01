/*!
   \file
   \brief Test class for cubic law permeability model.

   \copyright
    Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#include <gtest/gtest.h>

#include <memory>

#include "MaterialLib/MPL/Properties/PengRobinson.h"
#include "ParameterLib/ConstantParameter.h"
#include "TestMPL.h"

struct PengRobinsonTestData
{
    double pressure;          // p in Pa
    double temperature;       // T in K
    double expected_rho;      // rho in mol/m^3
    double expected_drho_dT;  // d(rho)/dT
    double expected_drho_dp;  // d(rho)/dp
};

// Table with test data points
const std::vector<PengRobinsonTestData> test_data = {
    {100000, 290, 1.83618, -0.00644593, 1.84721e-05},
    {120000, 290, 2.20607, -0.00777226, 1.85169e-05},
    {140000, 290, 2.57685, -0.00911135, 1.85619e-05},
    {160000, 290, 2.94854, -0.0104632, 1.8607e-05},
    {180000, 290, 3.32113, -0.0118284, 1.86525e-05},
    {200000, 290, 3.69464, -0.0132065, 1.8698e-05},
    {300000, 290, 5.57601, -0.0202995, 1.89305e-05},
    {500000, 290, 9.41016, -0.0355579, 1.94157e-05},
    {1e+06, 290, 19.4467, -0.0810047, 2.07654e-05},
    {100000, 300, 1.77393, -0.00600957, 1.78356e-05},
    {120000, 300, 2.13103, -0.00724441, 1.78741e-05},
    {140000, 300, 2.48891, -0.00848894, 1.79137e-05},
    {160000, 300, 2.84756, -0.0097445, 1.79526e-05},
    {180000, 300, 3.20701, -0.0110106, 1.7992e-05},
    {200000, 300, 3.56725, -0.0122888, 1.80321e-05},
    {300000, 300, 5.38042, -0.018847, 1.82324e-05},
    {500000, 300, 9.06836, -0.0328593, 1.86506e-05},
    {1e+06, 300, 18.6732, -0.0738831, 1.97951e-05},
    {100000, 330, 1.61042, -0.00494262, 1.6169e-05},
    {120000, 330, 1.93406, -0.00595108, 1.61948e-05},
    {140000, 330, 2.25822, -0.00696653, 1.6221e-05},
    {160000, 330, 2.5829, -0.0079886, 1.62473e-05},
    {180000, 330, 2.90811, -0.0090175, 1.62734e-05},
    {200000, 330, 3.23384, -0.0100531, 1.62998e-05},
    {300000, 330, 4.87047, -0.0153373, 1.64329e-05},
    {500000, 330, 8.18421, -0.0264507, 1.67059e-05},
    {1e+06, 330, 16.7155, -0.0577228, 1.74298e-05},
    {100000, 340, 1.56248, -0.00464999, 1.56819e-05},
    {120000, 340, 1.87635, -0.00559723, 1.57047e-05},
    {140000, 340, 2.19067, -0.00655027, 1.57275e-05},
    {160000, 340, 2.50545, -0.00750936, 1.57507e-05},
    {180000, 340, 2.82069, -0.00847416, 1.57737e-05},
    {200000, 340, 3.1364, -0.00944482, 1.57967e-05},
    {300000, 340, 4.72191, -0.0143892, 1.59136e-05},
    {500000, 340, 7.92838, -0.0247449, 1.6152e-05},
    {1e+06, 340, 16.1594, -0.0535893, 1.67796e-05},
    {100000, 360, 1.47475, -0.00413809, 1.47921e-05},
    {120000, 360, 1.77077, -0.00497892, 1.48097e-05},
    {140000, 360, 2.06714, -0.00582381, 1.48276e-05},
    {160000, 360, 2.36388, -0.00667333, 1.48456e-05},
    {180000, 360, 2.66097, -0.00752721, 1.48636e-05},
    {200000, 360, 2.95842, -0.00838551, 1.48815e-05},
    {300000, 360, 4.45108, -0.0127452, 1.4972e-05},
    {500000, 360, 7.46379, -0.0218114, 1.51558e-05},
    {1e+06, 360, 15.1598, -0.0466343, 1.56324e-05},
    {100000, 400, 1.32605, -0.00334035, 1.3288e-05},
    {120000, 400, 1.59192, -0.00401599, 1.32991e-05},
    {140000, 400, 1.85801, -0.00469439, 1.33101e-05},
    {160000, 400, 2.12432, -0.00537543, 1.33212e-05},
    {180000, 400, 2.39086, -0.00605903, 1.33323e-05},
    {200000, 400, 2.65761, -0.00674518, 1.33433e-05},
    {300000, 400, 3.99472, -0.0102158, 1.33988e-05},
    {500000, 400, 6.68568, -0.0173576, 1.35109e-05},
    {1e+06, 400, 13.512, -0.0364283, 1.37953e-05}};

TEST(MaterialPropertyLib, checkPengRobinsonEOS)
{
    const double Tc = 304.12;           // critical temperature (CO2) in K
    const double pc = 7.377e6;          // critical pressure (CO2) in Pa
    const double omega = 0.2249;        // acentric factor (CO2)
    const double molar_mass = 0.04401;  // molar_mass (CO2)

    MPL::Property const& density_model = MPL::PengRobinson{Tc, pc, omega};

    for (const auto& data : test_data)
    {
        // Create a VariableArray and set p, T, and molar_mass
        MPL::VariableArray variable_array;
        variable_array.gas_phase_pressure = data.pressure;
        variable_array.temperature = data.temperature;
        variable_array.molar_mass = molar_mass;

        ParameterLib::SpatialPosition const pos;
        double const t = std::numeric_limits<double>::quiet_NaN();
        double const dt = std::numeric_limits<double>::quiet_NaN();

        // Compute the density rho(p,T) using the Peng-Robinson model
        double rho =
            std::get<double>(density_model.value(variable_array, pos, t, dt));

        // Check if the calculated d(rho)/dT matches the expected value
        EXPECT_NEAR(data.expected_rho, rho, 1.e-4)
            << "Mismatch for rho at (p, T) = (" << data.pressure << ", "
            << data.temperature << ")";

        // Compute the derivative d(rho)/dT using the Peng-Robinson model
        double drho_dT = std::get<double>(density_model.dValue(
            variable_array, MPL::Variable::temperature, pos, t, dt));

        // Check if the calculated d(rho)/dT matches the expected value
        EXPECT_NEAR(data.expected_drho_dT, drho_dT, 1.e-4)
            << "Mismatch for d(rho)/dT at (p, T) = (" << data.pressure << ", "
            << data.temperature << ")";

        // Compute the derivative d(rho)/dp using the Peng-Robinson model
        double drho_dp = std::get<double>(density_model.dValue(
            variable_array, MPL::Variable::gas_phase_pressure, pos, t, dt));

        // Check if the calculated d(rho)/dp matches the expected value
        EXPECT_NEAR(data.expected_drho_dp, drho_dp, 1.e-4)
            << "Mismatch for d(rho)/dp at (p, T) = (" << data.pressure << ", "
            << data.temperature << ")";
    }
}
