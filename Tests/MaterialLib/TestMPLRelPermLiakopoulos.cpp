/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <gtest/gtest.h>

#include <sstream>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/RelativePermeability/RelPermLiakopoulos.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

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

    std::array<double, 17> ref_dk_rel_L_ds_L = {
        2.2276817495E+00, 2.2276817495E+00, 2.2276817495E+00, 2.2276817495E+00,
        2.2259427960E+00, 2.2240853238E+00, 2.2199407746E+00, 2.2150487817E+00,
        2.2090761349E+00, 2.2013998128E+00, 2.1965486778E+00, 2.1906259146E+00,
        2.1830137051E+00, 2.1723297907E+00, 2.1541864761E+00, 0.0000000000E+00,
        0.0000000000E+00};

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
