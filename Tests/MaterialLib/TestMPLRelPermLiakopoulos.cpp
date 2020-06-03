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

#include "TestMPL.h"
#include "Tests/TestTools.h"

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/RelativePermeability/RelPermLiakopoulos.h"

TEST(MaterialPropertyLib, RelPermLiakopoulos)
{
    std::stringstream m_beg;
    std::stringstream m_end;

    m_beg << "<medium>\n";
    m_beg << "<phases></phases>\n";
    m_beg << "<properties>\n";
    m_beg << "  <property>\n";
    m_beg << "    <name>saturation</name>\n";
    m_beg << "    <type>Constant</type>\n";
    // constant saturation value will be inserted here...
    m_end << "  </property>\n";
    m_end << "  <property>\n";
    m_end << "    <name>relative_permeability</name>\n";
    m_end << "    <type>RelPermLiakopoulos</type>\n";
    m_end << "  </property>\n";
    m_end << "</properties>\n";
    m_end << "</medium>\n";

    std::array<double, 17> ref_saturation = {0.00, 0.10, 0.15, 0.20, 0.25, 0.30,
                                             0.40, 0.50, 0.60, 0.70, 0.75, 0.80,
                                             0.85, 0.90, 0.95, 1.00, 1.05};

    std::array<double, 17> ref_k_rel_L = {
        0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00,
        0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00,
        1.2693364888E-01, 3.4747560139E-01, 4.5742795233E-01, 5.6711275277E-01,
        6.7646274501E-01, 7.8536411514E-01, 8.9357837782E-01, 1.0000000000E+00,
        1.0000000000E+00};

    std::array<double, 17> ref_k_rel_G = {
        1.0000000000E+00, 1.0000000000E+00, 1.0000000000E+00, 1.0000000000E+00,
        8.7025505871E-01, 7.4169921875E-01, 5.0669293177E-01, 3.1445005531E-01,
        1.7125493438E-01, 7.6376563369E-02, 4.5357912805E-02, 2.3805540052E-02,
        1.0284395970E-02, 3.1176074364E-03, 3.9836380435E-04, 1.0000000000E-04,
        1.0000000000E-04};

    std::array<double, 17> ref_dk_rel_L_ds_L = {
        2.2276817495E+00, 2.2276817495E+00, 2.2276817495E+00, 2.2276817495E+00,
        2.2259427960E+00, 2.2240853238E+00, 2.2199407746E+00, 2.2150487817E+00,
        2.2090761349E+00, 2.2013998128E+00, 2.1965486778E+00, 2.1906259146E+00,
        2.1830137051E+00, 2.1723297907E+00, 2.1541864761E+00, 0.0000000000E+00,
        0.0000000000E+00};

    std::array<double, 17> ref_dk_rel_G_ds_L = {
        -2.5000000000E+00, -2.5000000000E+00, -2.5000000000E+00,
        -2.5000000000E+00, -2.6090531996E+00, -2.5179036458E+00,
        -2.1540353412E+00, -1.6809943584E+00, -1.1843791120E+00,
        -7.2333854456E-01, -5.2134311212E-01, -3.4554001149E-01,
        -2.0089926430E-01, -9.2131654833E-02, -2.3729854831E-02,
        0.0000000000E+00,  0.0000000000E+00};

    for (size_t idx = 0; idx < ref_saturation.size(); idx++)
    {
        std::stringstream m_sat;
        m_sat << "    <value>" << ref_saturation[idx] << "</value>\n";
        std::stringstream m;
        m << m_beg.str() << m_sat.str() << m_end.str();

        auto const& medium = Tests::createTestMaterial(m.str());
        MaterialPropertyLib::VariableArray variable_array;
        ParameterLib::SpatialPosition const pos;
        double const time = std::numeric_limits<double>::quiet_NaN();
        double const dt = std::numeric_limits<double>::quiet_NaN();

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
