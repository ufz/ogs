/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationVanGenuchtenWithVolumetricStrain.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

namespace MPL = MaterialPropertyLib;

TEST(MaterialPropertyLib, SaturationVanGenuchtenWithVolumetricStrain)
{
    double const residual_liquid_saturation = 0.1;
    double const residual_gas_saturation = 0.05;
    double const exponent = 0.5;
    double const p_b = 5000000;
    double const e_0 = 0.7;
    double const e_m = 0.5;
    double const a = 8;
    double const d_diff = 20;

    MPL::Property const& pressure_saturation =
        MPL::SaturationVanGenuchtenWithVolumetricStrain{
            "saturation",
            residual_liquid_saturation,
            residual_gas_saturation,
            exponent,
            p_b,
            e_0,
            e_m,
            a,
            d_diff};

    MPL::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    /// No volumetric strain
    {
        double const vol_s = 0;
        double const p_L = 1000000;

        variable_array.capillary_pressure = p_L;
        variable_array.volumetric_strain = vol_s;

        double const S = pressure_saturation.template value<double>(
            variable_array, pos, t, dt);

        double const S_expected = 0.75426406;

        ASSERT_NEAR(S, S_expected, 1e-5);
    }
    /// Volumetric strain
    {
        double const vol_s = 0.002;
        double const p_L = 1000000;

        variable_array.capillary_pressure = p_L;
        variable_array.volumetric_strain = vol_s;

        double const S = pressure_saturation.template value<double>(
            variable_array, pos, t, dt);

        double const S_expected = 0.73068816;

        ASSERT_NEAR(S, S_expected, 1e-5);
    }
    // Test derivative (dValue) of saturation with respect to capillary pressure
    {
        double const vol_s = 0.002;
        double const p_L = 1000000;

        variable_array.capillary_pressure = p_L;
        variable_array.volumetric_strain = vol_s;

        // Calculate saturation
        double const S = pressure_saturation.template value<double>(
            variable_array, pos, t, dt);

        // Calculate the derivative numerically using a small perturbation
        double const dP = 1.0e-6;
        variable_array.capillary_pressure = p_L + dP;
        double const S_plus_dP = pressure_saturation.template value<double>(
            variable_array, pos, t, dt);

        double const approximated_dS_dP = (S_plus_dP - S) / dP;

        // Calculate the derivative analytically using dValue
        double const analytic_dS_dP =
            pressure_saturation.template dValue<double>(
                variable_array,
                MaterialPropertyLib::Variable::capillary_pressure, pos, t, dt);

        // Check if the numerical and analytical derivatives are close
        ASSERT_LE(std::fabs(approximated_dS_dP - analytic_dS_dP), 1e-6)
            << "for expected derivative of saturation with respect to "
               "capillary pressure "
            << approximated_dS_dP
            << " and for computed derivative of saturation with respect to "
               "capillary pressure "
            << analytic_dS_dP;
    }
}
