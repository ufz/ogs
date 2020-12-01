/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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

TEST(MaterialPropertyLib, RelPermUdell)
{
    const double ref_residual_liquid_saturation = 0.01;
    const double ref_residual_gas_saturation = 0.01;
    const double ref_k_rel_L_min = 1.e-9;
    const double ref_k_rel_G_min = 1.e-7;

    std::stringstream m;

    m << "<medium>\n";
    m << "<properties>\n";
    m << "  <property>\n";
    m << "    <name>relative_permeability</name>\n";
    m << "    <type>RelativePermeabilityUdell</type>\n";
    m << "    <residual_liquid_saturation>" << ref_residual_liquid_saturation
      << "</residual_liquid_saturation>\n";
    m << "    <residual_gas_saturation>" << ref_residual_gas_saturation
      << "</residual_gas_saturation>\n";
    m << "    <min_relative_permeability_liquid>" << ref_k_rel_L_min
      << "</k_rel_min_liquid>\n";
    m << "    <min_relative_permeability_gas>" << ref_k_rel_G_min
      << "</k_rel_min_gas>\n";
    m << "  </property>\n";
    m << "</properties>\n";
    m << "</medium>\n";

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

    std::array<double, 16> ref_k_rel_G = {
        1.00000000000E+00, 1.00000000000E+00, 1.00000000000E+00,
        9.69699062465E-01, 9.10945906893E-01, 7.49017203716E-01,
        5.23845294053E-01, 2.18211587009E-01, 6.30253975809E-02,
        7.28756725514E-03, 7.74549719930E-04, 2.86870266641E-05,
        1.06248246904E-06, 1.00000000000E-07, 1.00000000000E-07,
        1.00000000000E-07};

    std::array<double, 16> ref_dk_rel_L_ds_L = {
        0.00000000000E+00, 0.00000000000E+00, 0.00000000000E+00,
        3.18744740712E-04, 2.86870266641E-03, 2.58183239977E-02,
        1.15066851397E-01, 4.84810750623E-01, 1.10955044242E+00,
        1.98928592678E+00, 2.52477709118E+00, 2.87667128492E+00,
        2.99906926536E+00, 3.06122448980E+00, 0.00000000000E+00,
        0.00000000000E+00};

    std::array<double, 16> ref_dk_rel_G_ds_L = {
        0.00000000000E+00, 0.00000000000E+00, 3.06122448980E+00,
        2.99906926536E+00, 2.87667128492E+00, 2.52477709118E+00,
        1.98928592678E+00, 1.10955044242E+00, 4.84810750623E-01,
        1.15066851397E-01, 2.58183239977E-02, 2.86870266641E-03,
        3.18744740712E-04, 0.00000000000E+00, 0.00000000000E+00,
        0.00000000000E+00};

    for (std::size_t idx = 0; idx < ref_saturation.size(); idx++)
    {
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
                .template value<Eigen::Vector2d>(variable_array, pos, time, dt);

        auto dk_rel_ds_L =
            medium
                ->property(
                    MaterialPropertyLib::PropertyType::relative_permeability)
                .template dValue<Eigen::Vector2d>(
                    variable_array,
                    MaterialPropertyLib::Variable::liquid_saturation, pos, time,
                    dt);

        auto k_rel_L = k_rel[0];
        auto k_rel_G = k_rel[1];

        auto dk_rel_L_ds_L = dk_rel_ds_L[0];
        auto dk_rel_G_ds_L = dk_rel_ds_L[1];

        ASSERT_NEAR(k_rel_L, ref_k_rel_L[idx], 1.0e-10);
        ASSERT_NEAR(k_rel_G, ref_k_rel_G[idx], 1.0e-10);
        ASSERT_NEAR(dk_rel_L_ds_L, ref_dk_rel_L_ds_L[idx], 1.0e-10);
        ASSERT_NEAR(dk_rel_G_ds_L, ref_dk_rel_G_ds_L[idx], 1.0e-10);
    }
}
