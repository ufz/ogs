// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <map>
#include <range/v3/view/enumerate.hpp>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationLuMcCartney.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

namespace MPL = MaterialPropertyLib;

std::stringstream configSaturationLuMcCartney(std::string const& material)
{
    std::string material_name = material;
    std::stringstream m;
    double rho_s;
    double phi;
    if (material == "BoomClay")
    {
        phi = 0.483;
        rho_s = 1400.0 / (1 - phi);
    }
    else if (material == "FEBEX")
    {
        phi = 0.443;
        rho_s = 1500.0 / (1 - phi);
    }
    else if (material == "MX80")
    {
        phi = 0.432;
        rho_s = 1600 / (1 - phi);
    }
    else if (material == "GMZ01")
    {
        phi = 0.359;
        rho_s = 1700.0 / (1 - phi);
    }
    else
    {
        material_name = "MX80";
        phi = 0.331;
        rho_s = 2500.0;
    }

    m << "<medium>\n";
    m << "<phases>\n";
    m << "  <phase>\n";
    m << "    <type>Solid</type>\n";
    m << "    <properties>\n";
    m << "        <property>\n";
    m << "            <name>density</name>\n";
    m << "            <type>Constant</type>\n";
    m << "            <value>" << rho_s << "</value>\n";
    m << "        </property>\n";
    m << "    </properties>\n";
    m << "  </phase>\n";
    m << "</phases>\n";
    m << "<properties>\n";
    m << "    <property>\n";
    m << "        <name>porosity</name>\n";
    m << "        <type>Constant</type>\n";
    m << "        <value>" << phi << "</value>\n";
    m << "    </property>\n";
    m << "    <property>\n";
    m << "        <name>saturation</name>\n";
    m << "        <type>SaturationLuMcCartney</type>\n";
    m << "        <material>" << material_name << "</material>\n";
    m << "    </property>\n";
    m << "</properties>\n";
    m << "</medium>\n";
    return m;
}

TEST(MaterialPropertyLib, SaturationLuMcCartney)
{
    std::array<double, 3> test_pressures = {113.0e6, 38.0e6, 4200e3};
    std::map<std::string, std::array<double, 2>> test_temperatures;
    test_temperatures["BoomClay"] = {22, 80};
    test_temperatures["FEBEX"] = {26, 80};
    test_temperatures["MX80"] = {20, 60};
    test_temperatures["GMZ01"] = {20, 80};
    std::map<std::string, std::array<std::array<double, 2>, 3>>
        ref_saturation_values;
    ref_saturation_values["BoomClay"] = {
        {{0.07729657544222532, 0.059762990733618036},
         {0.16424492073037383, 0.1452509769618063},
         {0.3703627576499768, 0.3700666607141584}}};
    ref_saturation_values["FEBEX"] = {
        {{0.4591891631738408, 0.37699611659111154},
         {0.6527266897490632, 0.6337364514817921},
         {0.8389461121833218, 0.8267385116139065}}};
    ref_saturation_values["MX80"] = {
        {{0.34391268352030857, 0.2878185477378836},
         {0.5717438358984998, 0.5434054756227413},
         {0.9261226316345796, 0.9110943287569135}}};
    ref_saturation_values["GMZ01"] = {
        {{0.5684361050669597, 0.45948154104341293},
         {0.720296958963634, 0.7048034727066785},
         {0.9618659666552256, 0.9144921134069371}}};
    std::map<std::string, std::stringstream> m;
    m["MX80_Bentonite"] = configSaturationLuMcCartney("MX80_Bentonite");

    for (auto const& entry : test_temperatures)
    {
        m[entry.first] = configSaturationLuMcCartney(entry.first);
    }
    auto medium = Tests::createTestMaterial(m["MX80"].str());

    MPL::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    auto S =
        [pos, t, dt, &medium, &variable_array](double const p_L, double const T)
    {
        variable_array.capillary_pressure = -p_L;
        variable_array.temperature = T;
        return medium->property(MaterialPropertyLib::PropertyType::saturation)
            .template value<double>(variable_array, pos, t, dt);
    };

    // Test reference values
    for (auto const& entry : test_temperatures)
    {
        for (auto const& [i, p] : ranges::views::enumerate(test_pressures))
        {
            for (auto const& [j, T] :
                 ranges::views::enumerate(test_temperatures[entry.first]))
            {
                medium = Tests::createTestMaterial(m[entry.first].str());
                double const S_L = S(-p, T + 273.15);
                // The reference values are regenerated at full double
                // precision from the current constants, so this band buys
                // only numerically slack: it absorbs the realisation noise
                // with which another maths library runs the model's
                // exp/pow/erf chain - a few ulps of S per evaluation
                // (2.1e-17 to 3.4e-16 at these saturations), amplified
                // along the chain. The 5e-12 band stands more than 1e4
                // times the few-ulp ceiling and admits real drift up to
                // 5e-12; drift beyond the band reads red.
                EXPECT_LE(
                    std::abs(S_L - ref_saturation_values[entry.first][i][j]),
                    5e-12)
                    << "for temperature " << T << ", capillary pressure " << p
                    << " and saturation " << S(-p, T + 273.15) << "\n";
            }
        }
    }

    // Test derivatives
    double const p_0 = -1e9;
    double const p_max = 1e6;
    double const T_0 = 273.15;
    double const T_max = 120.0 + T_0;
    int const p_steps = 10000;
    int const T_steps = 10;
    // Fraction of total steps allocated to positive pressure range [0, p_max]
    double const positive_fraction = 0.1;
    // Concentration factor: higher values = more points near p_L = 0
    double const alpha = 2.0;
    medium = Tests::createTestMaterial(m["MX80_Bentonite"].str());

    for (int i = 0; i <= p_steps; ++i)
    {
        double const lambda = static_cast<double>(i) / p_steps;  // [0, 1]
        for (int j = 0; j <= T_steps; ++j)
        {
            double p_L;
            if (lambda <= 1.0 - positive_fraction)
            {
                // Negative pressure region: [p_0, 0]
                // Power function concentrates points near p_L = 0
                double const s = lambda / (1.0 - positive_fraction);  // [0, 1]
                p_L = p_0 * std::pow(1.0 - s, alpha);
            }
            else
            {
                // Positive pressure region: [0, p_max]
                double const s = (lambda - (1.0 - positive_fraction)) /
                                 positive_fraction;  // [0, 1]
                p_L = p_max * s;
            }

            double const T = T_0 + j * (T_max - T_0) / T_steps;

            double const S_L = S(p_L, T);
            double const dS_p =
                medium->property(MaterialPropertyLib::PropertyType::saturation)
                    .template dValue<double>(variable_array,
                                             MPL::Variable::capillary_pressure,
                                             pos, t, dt);

            double const dS2_p =
                medium->property(MaterialPropertyLib::PropertyType::saturation)
                    .template d2Value<double>(
                        variable_array, MPL::Variable::capillary_pressure,
                        MPL::Variable::capillary_pressure, pos, t, dt);

            double const dS_T =
                medium->property(MaterialPropertyLib::PropertyType::saturation)
                    .template dValue<double>(
                        variable_array, MPL::Variable::temperature, pos, t, dt);

            double const dS2_T =
                medium->property(MaterialPropertyLib::PropertyType::saturation)
                    .template d2Value<double>(
                        variable_array, MPL::Variable::temperature,
                        MPL::Variable::temperature, pos, t, dt);
            double const dS2_pT =
                medium->property(MaterialPropertyLib::PropertyType::saturation)
                    .template d2Value<double>(
                        variable_array, MPL::Variable::capillary_pressure,
                        MPL::Variable::temperature, pos, t, dt);

            // Each comparison band below covers the worst-case round-off
            // of its own finite-difference stencil: a few ulps of S
            // (ulp(S) <= ulp(1) = 2.2e-16 at these saturations) per
            // evaluation, amplified by the stencil's coefficient weights
            // over its denominator, with no cancellation between the
            // evaluations assumed (the second-difference (1, -2, 1)
            // weights amplify aligned error instead). The tightest case is
            // the 5e-7 dS2_T band, where a three-ulp realisation leaves
            // only about a factor of two; the remaining bands sit roughly
            // four times or more above their own worst case.
            double eps = 1.0e-2;
            double const DS_p = (S_L - S(p_L + eps, T)) / (eps);
            EXPECT_LE(std::abs(dS_p - DS_p), 5e-13)
                << "for temperature " << T << ", capillary pressure " << -p_L
                << " and saturation " << S_L << " dS_p " << dS_p << " DS_p "
                << DS_p << " eps " << eps << "\n";

            eps = 100;
            double const DS2_p =
                (S(p_L + eps, T) - 2 * S(p_L, T) + S(p_L - eps, T)) /
                (eps * eps);

            EXPECT_LE(std::abs(dS2_p - DS2_p), 1e-10)
                << "for temperature " << T << ", capillary pressure " << -p_L
                << " and saturation " << S_L << " dS2_p " << dS2_p << " DS2_p "
                << DS2_p << " eps " << eps << "\n";

            eps = 1e-4;
            double const DS_T = (S(p_L, T + eps) - S(p_L, T - eps)) / (2 * eps);
            // double const DS_T = (S(p_L, T) - S(p_L, T - eps)) / (eps);
            EXPECT_LE(std::abs(dS_T - DS_T), 5e-11)
                << "for temperature " << T << ", capillary pressure " << -p_L
                << " and saturation " << S_L << " dS_T " << dS_T << " DS_T "
                << DS_T << " eps " << eps << "\n";

            eps = 1e-4;
            double const DS2_T =
                (S(p_L, T + eps) - 2 * S(p_L, T) + S(p_L, T - eps)) /
                (eps * eps);
            EXPECT_LE(std::abs(dS2_T - DS2_T), 5e-7)
                << "for temperature " << T << ", capillary pressure " << -p_L
                << " and saturation " << S_L << " dS2_T " << dS2_T << " DS2_T "
                << DS2_T << " eps " << eps << "\n";

            double const eps_p = 1;
            double const DS2_pT =
                (S(p_L + eps_p, T - eps) + S(p_L - eps_p, T + eps) -
                 S(p_L - eps_p, T - eps) - S(p_L + eps_p, T + eps)) /
                (4 * eps_p * eps);
            EXPECT_LE(std::abs(dS2_pT - DS2_pT), 5e-11)
                << "for temperature " << T << ", capillary pressure " << -p_L
                << " and saturation " << S(p_L, T) << " dS2_pT " << dS2_pT
                << " DS2_pT " << DS2_pT << " eps " << eps << " eps_p " << eps_p
                << "\n";
        }
    }
}
