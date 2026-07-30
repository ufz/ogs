// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <limits>

#include "MaterialLib/MPL/Properties/Density/CreateWaterVapourDensity.h"
#include "TestMPL.h"

static double centralDifferencesDerivative(
    double const value,
    double MaterialPropertyLib::VariableArray::* const variable,
    double const value_increment, MaterialPropertyLib::Property const& property,
    MaterialPropertyLib::VariableArray variable_array, const auto& pos,
    double const t, double const dt)
{
    variable_array.*variable = value - value_increment;
    double const value_minus =
        property.template value<double>(variable_array, pos, t, dt);

    variable_array.*variable = value + value_increment;
    double const value_plus =
        property.template value<double>(variable_array, pos, t, dt);

    return 0.5 * (value_plus - value_minus) / value_increment;
}

TEST(MaterialPropertyLib, WaterVapourDensity)
{
    char const xml[] =
        "<property>"
        "   <name>vapour_density</name>"
        "   <type>WaterVapourDensity</type>"
        "</property>";

    std::unique_ptr<MaterialPropertyLib::Property> const property_ptr =
        Tests::createTestProperty(
            xml, MaterialPropertyLib::createWaterVapourDensity);
    MaterialPropertyLib::Property const& property = *property_ptr;

    double const T = 390.0;
    double const p = -1.e+6;
    double const rho_w = 1000.0;

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    variable_array.temperature = T;
    variable_array.liquid_phase_pressure = p;
    variable_array.density = rho_w;

    // The derivative of the water vapour with respect of temperature
    {
        // The first evaluation point 273.0001 is chosen to avoid lower limit
        // (T=273K) of the function in when computing numerical derivative.
        std::array const temperatures = {273.0001, 293.0, 393.0, 420.0, 500.0};
        std::array const rho_vw_expected = {
            4.8760232545732324e-03, 1.6928712260750413e-02, 1.2768654862482858,
            2.8826351745119623, 1.9203858907405188e+01};

        for (std::size_t i = 0; i < temperatures.size(); ++i)
        {
            double const T_i = temperatures[i];
            variable_array.temperature = T_i;

            double const rho_vw =
                property.template value<double>(variable_array, pos, t, dt);

            // The expected values are re-derived from the documented
            // formula at full double precision, so this band buys only
            // numerically slack: it absorbs the realisation noise with
            // which another maths library runs the correlation's exp
            // chain - a few ulps of the result per evaluation (at most
            // 1.1e-14 at the largest result, 19.2), amplified along the
            // chain. The 5e-12 band stands more than 400 times that
            // per-evaluation ceiling and admits real drift up to 5e-12;
            // drift beyond the band reads red.
            ASSERT_LE(std::fabs(rho_vw_expected[i] - rho_vw), 5e-12)
                << "for expected water vapour density " << rho_vw_expected[i]
                << " and for computed water vapour density " << rho_vw;

            double const approximated_drho_wv_dT = centralDifferencesDerivative(
                T_i, &MaterialPropertyLib::VariableArray::temperature, 1e-4,
                property, variable_array, pos, t, dt);

            double const analytic_drho_wv_dT = property.template dValue<double>(
                variable_array, MaterialPropertyLib::Variable::temperature, pos,
                t, dt);

            // The increment of 1e-4 K is the largest one keeping T - eps at
            // or above the 273 K evaluation limit at the first evaluation
            // point. The observed analytic-vs-stencil gap reaches 1.0e-10 at
            // T = 500 K, so the band is 5e-10: about a factor of five above
            // the observed gap, with room for the few ulps per evaluation
            // that another realisation of std::exp may introduce.
            EXPECT_LE(std::fabs(approximated_drho_wv_dT - analytic_drho_wv_dT),
                      5e-10)
                << "for expected derivative of water vapour density with "
                   "respect to temperature "
                << approximated_drho_wv_dT
                << " and for computed derivative of water vapour density with "
                   "respect to temperature."
                << analytic_drho_wv_dT;
        }
    }

    // The derivative of the water vapour density with respect of pressure
    {
        std::array const pressures = {-1.e+7, -1.e+6, 1.e+5, 1.e+6,
                                      2.e+6,  6.e+6,  1.e+7};
        std::array const rho_vw_expected = {
            1.1018266236068592, 1.1583204083416694, 1.1654209431266414,
            1.1712628346385106, 1.177788176514993,  1.2042551150024834,
            1.2313168113987973};

        for (std::size_t i = 0; i < pressures.size(); ++i)
        {
            double const p_i = pressures[i];

            variable_array.temperature = T;
            variable_array.liquid_phase_pressure = p_i;

            double const rho_vw =
                property.template value<double>(variable_array, pos, t, dt);

            // Re-derived at full double precision as above; a few ulps
            // of the result is below 7e-16 here, and drift beyond the
            // 5e-12 band reads red.
            ASSERT_LE(std::fabs(rho_vw_expected[i] - rho_vw), 5e-12)
                << "for expected water vapour density " << rho_vw_expected[i]
                << " and for computed water vapour density " << rho_vw;

            // A 1 Pa increment keeps the stencil's rounding error at the
            // 1e-16 level: the correlation chains two std::exp calls with
            // arithmetic and can err by a few ulps of the result per
            // evaluation, so the 2 Pa denominator admits a worst case near
            // 1e-15 at these magnitudes - with no cancellation between the
            // two evaluations assumed - while the truncation error stays
            // negligible: rho(p) is linear in p to well below 1e-15
            // relative over the +/-1 Pa stencil. The 2e-13 band keeps more
            // than two orders of magnitude of margin either way.
            double const approximated_drho_wv_dp = centralDifferencesDerivative(
                p_i, &MaterialPropertyLib::VariableArray::liquid_phase_pressure,
                1.0, property, variable_array, pos, t, dt);

            double const analytic_drho_wv_dp = property.template dValue<double>(
                variable_array,
                MaterialPropertyLib::Variable::liquid_phase_pressure, pos, t,
                dt);

            EXPECT_LE(std::fabs(approximated_drho_wv_dp - analytic_drho_wv_dp),
                      2e-13)
                << "for expected derivative of water vapour density with "
                   "respect to pressure "
                << approximated_drho_wv_dp
                << " and for computed derivative of water vapour density "
                   "with respect to pressure."
                << analytic_drho_wv_dp;
        }
    }
}
