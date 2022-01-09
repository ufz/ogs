/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <gtest/gtest.h>

#include <sstream>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/RelativePermeability/RelPermUdell.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

std::stringstream configRelPermUdell(
    const double ref_residual_liquid_saturation,
    const double ref_residual_gas_saturation,
    const double ref_k_rel_min)
{
    std::stringstream m;

    m << "<medium>\n";
    m << "<properties>\n";
    m << "  <property>\n";
    m << "    <name>relative_permeability</name>\n";
    m << "    <type>RelativePermeabilityUdell</type>\n";
    m << "    <residual_liquid_saturation>" << ref_residual_liquid_saturation
      << "    </residual_liquid_saturation>\n";
    m << "    <residual_gas_saturation>" << ref_residual_gas_saturation
      << "    </residual_gas_saturation>\n";
    m << "    <min_relative_permeability>" << ref_k_rel_min
      << "    </min_relative_permeability>\n";
    m << "  </property>\n";
    m << "</properties>\n";
    m << "</medium>\n";
    return m;
}

TEST(MaterialPropertyLib, RelPermUdellDerivatives)
{
    const double ref_residual_liquid_saturation = 0.01;
    const double ref_residual_gas_saturation = 0.01;
    const double ref_k_rel_L_min = 1.e-9;
    const double s_max = 1. - ref_residual_gas_saturation;
    const double stepsize = 0.01;

    for (double s_L = ref_residual_liquid_saturation + stepsize;
         s_L <= s_max - stepsize;
         s_L += stepsize)
    {
        // Wetting phase
        const std::stringstream m =
            configRelPermUdell(ref_residual_liquid_saturation,
                               ref_residual_gas_saturation, ref_k_rel_L_min);

        auto const& medium = Tests::createTestMaterial(m.str());

        MaterialPropertyLib::VariableArray variable_array;
        ParameterLib::SpatialPosition const pos;
        double const time = std::numeric_limits<double>::quiet_NaN();
        double const dt = std::numeric_limits<double>::quiet_NaN();

        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::liquid_saturation)] = s_L;

        auto const k_rel =
            medium
                ->property(
                    MaterialPropertyLib::PropertyType::relative_permeability)
                .template value<double>(variable_array, pos, time, dt);

        auto const dk_rel_ds_L =
            medium
                ->property(
                    MaterialPropertyLib::PropertyType::relative_permeability)
                .template dValue<double>(
                    variable_array,
                    MaterialPropertyLib::Variable::liquid_saturation, pos, time,
                    dt);

        const double eps = 1.e-8;

        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::liquid_saturation)] = s_L + eps;

        auto const k_rel_plus =
            medium
                ->property(
                    MaterialPropertyLib::PropertyType::relative_permeability)
                .template value<double>(variable_array, pos, time, dt);

        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::liquid_saturation)] = s_L - eps;
        auto const k_rel_minus =
            medium
                ->property(
                    MaterialPropertyLib::PropertyType::relative_permeability)
                .template value<double>(variable_array, pos, time, dt);

        const auto Dk_rel = (k_rel_plus - k_rel_minus) / (2 * eps);

        if ((s_L < ref_residual_liquid_saturation) || (s_L > s_max))
        {
            ASSERT_EQ(Dk_rel, 0.)
                << "for saturation " << s_L
                << " and relative liquid phase permeability " << k_rel;
        }
        else
        {
            ASSERT_LE(std::abs(dk_rel_ds_L - Dk_rel), 1e-7)
                << "with " << dk_rel_ds_L << " and " << Dk_rel
                << " for saturation " << s_L
                << " and relative liquid phase permeability " << k_rel;
        }
    }
}

TEST(MaterialPropertyLib, RelPermUdellRefValues)
{
    std::array<double, 16> ref_saturation = {
        -0.01, 0.00, 0.01, 0.02, 0.04, 0.10, 0.20, 0.40,
        0.60,  0.80, 0.90, 0.96, 0.98, 0.99, 1.00, 1.01};

    std::array<double, 16> ref_k_rel_L = {
        1.00000000000E-09, 1.00000000000E-09, 1.00000000000E-09,
        1.06248246904E-06, 2.86870266641E-05, 7.74549719930E-04,
        7.28756725514E-03, 6.30253975809E-02, 2.18211587009E-01,
        5.23845294053E-01, 7.49017203716E-01, 9.10945906893E-01,
        9.69699062465E-01, 1.00000000000E+00, 1.00000000000E+00,
        1.00000000000E+00};

    std::array<double, 16> ref_dk_rel_L_ds_L = {
        0.00000000000E+00, 0.00000000000E+00, 0.00000000000E+00,
        3.18744740712E-04, 2.86870266641E-03, 2.58183239977E-02,
        1.15066851397E-01, 4.84810750623E-01, 1.10955044242E+00,
        1.98928592678E+00, 2.52477709118E+00, 2.87667128492E+00,
        2.99906926536E+00, 3.06122448980E+00, 0.00000000000E+00,
        0.00000000000E+00};

    for (std::size_t idx = 0; idx < ref_saturation.size(); idx++)
    {
        const double ref_residual_liquid_saturation = 0.01;
        const double ref_residual_gas_saturation = 0.01;
        const double ref_k_rel_L_min = 1.e-9;

        const std::stringstream m =
            configRelPermUdell(ref_residual_liquid_saturation,
                               ref_residual_gas_saturation, ref_k_rel_L_min);

        auto const& medium = Tests::createTestMaterial(m.str());
        MaterialPropertyLib::VariableArray variable_array;
        ParameterLib::SpatialPosition const pos;
        double const time = std::numeric_limits<double>::quiet_NaN();
        double const dt = std::numeric_limits<double>::quiet_NaN();

        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::liquid_saturation)] =
            ref_saturation[idx];

        auto k_rel =
            medium
                ->property(
                    MaterialPropertyLib::PropertyType::relative_permeability)
                .template value<double>(variable_array, pos, time, dt);

        auto dk_rel_ds_L =
            medium
                ->property(
                    MaterialPropertyLib::PropertyType::relative_permeability)
                .template dValue<double>(
                    variable_array,
                    MaterialPropertyLib::Variable::liquid_saturation, pos, time,
                    dt);

        ASSERT_NEAR(k_rel, ref_k_rel_L[idx], 1.0e-10);
        ASSERT_NEAR(dk_rel_ds_L, ref_dk_rel_L_ds_L[idx], 1.0e-10);
    }
}
