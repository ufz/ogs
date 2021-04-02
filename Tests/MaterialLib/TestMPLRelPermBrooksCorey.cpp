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
#include "MaterialLib/MPL/Properties/RelativePermeability/RelPermBrooksCorey.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

TEST(MaterialPropertyLib, RelPermBrooksCorey)
{
    const double ref_lambda = 2.5;
    const double ref_residual_liquid_saturation = 0.01;
    const double ref_residual_gas_saturation = 0.01;
    const double ref_k_rel_L_min = 1.e-9;

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
    m_end << "    <type>RelPermBrooksCorey</type>\n";
    m_end << "    <residual_liquid_saturation>"
          << ref_residual_liquid_saturation
          << "    </residual_liquid_saturation>\n";
    m_end << "    <residual_gas_saturation>" << ref_residual_gas_saturation
          << "    </residual_gas_saturation>\n";
    m_end << "    <lambda>" << ref_lambda << "</lambda>\n";
    m_end << "    <min_relative_permeability>" << ref_k_rel_L_min
          << "    </min_relative_permeability>\n";
    m_end << "  </property>\n";
    m_end << "</properties>\n";
    m_end << "</medium>\n";

    std::array<double, 16> ref_saturation = {
        -0.01, 0.00, 0.01, 0.02, 0.04, 0.10, 0.20, 0.40,
        0.60,  0.80, 0.90, 0.96, 0.98, 0.99, 1.00, 1.01};

    std::array<double, 16> ref_k_rel_L = {
        1.0000000000E-09, 1.0000000000E-09, 1.0000000000E-09, 2.7123199126E-08,
        1.7636064574E-06, 1.1467333635E-04, 1.9615735463E-03, 3.0156880778E-02,
        1.4540473747E-01, 4.4088351375E-01, 6.9346243084E-01, 8.8856788446E-01,
        9.6177504114E-01, 1.0000000000E+00, 1.0000000000E+00, 1.0000000000E+00};

    std::array<double, 16> ref_dk_rel_L_ds_L = {
        0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 1.0306815668E-05,
        2.2339015126E-04, 4.8417630905E-03, 3.9231470925E-02, 2.9383627425E-01,
        9.3650508877E-01, 2.1207055092E+00, 2.9608508283E+00, 3.5542715378E+00,
        3.7677785117E+00, 3.8775510204E+00, 0.0000000000E+00, 0.0000000000E+00};

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
