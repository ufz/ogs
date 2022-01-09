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
#include "MaterialLib/MPL/Properties/RelativePermeability/RelPermBrooksCoreyNonwettingPhase.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

TEST(MaterialPropertyLib, RelPermBrooksCoreyNonwettingPhase)
{
    const double ref_lambda = 2.5;
    const double ref_residual_liquid_saturation = 0.01;
    const double ref_residual_gas_saturation = 0.01;
    const double ref_k_rel_G_min = 1.e-7;

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
    m_end << "    <name>relative_permeability_nonwetting_phase</name>\n";
    m_end << "    <type>RelPermBrooksCoreyNonwettingPhase</type>\n";
    m_end << "    <residual_liquid_saturation>"
          << ref_residual_liquid_saturation
          << "    </residual_liquid_saturation>\n";
    m_end << "    <residual_gas_saturation>" << ref_residual_gas_saturation
          << "    </residual_gas_saturation>\n";
    m_end << "    <lambda>" << ref_lambda << "</lambda>\n";
    m_end << "    <min_relative_permeability>" << ref_k_rel_G_min
          << "    </min_relative_permeability>\n";
    m_end << "  </property>\n";
    m_end << "</properties>\n";
    m_end << "</medium>\n";

    std::array<double, 16> ref_saturation = {
        -0.01, 0.00, 0.01, 0.02, 0.04, 0.10, 0.20, 0.40,
        0.60,  0.80, 0.90, 0.96, 0.98, 0.99, 1.00, 1.01};

    std::array<double, 16> ref_k_rel_G = {
        1.0000000000E+00, 1.0000000000E+00, 1.0000000000E+00, 9.7944075784E-01,
        9.3794411438E-01, 8.1354659673E-01, 6.1592154540E-01, 2.9343532599E-01,
        9.4837870477E-02, 1.2086349932E-02, 1.3426518034E-03, 5.1003060118E-05,
        1.9046571241E-06, 1.0000000000E-07, 1.0000000000E-07, 1.0000000000E-07};

    std::array<double, 16> ref_dk_rel_G_ds_L = {
        0.0000000000E+00,  0.0000000000E+00,  -2.0408163265E+00,
        -2.0654018726E+00, -2.0807295100E+00, -2.0524729938E+00,
        -1.8805652791E+00, -1.3132397982E+00, -6.8017950205E-01,
        -1.8533091175E-01, -4.4178730635E-02, -5.0791425970E-03,
        -5.7061547092E-04, -0.0000000000E+00, 0.0000000000E+00,
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

        auto k_rel = medium
                         ->property(MaterialPropertyLib::PropertyType::
                                        relative_permeability_nonwetting_phase)
                         .template value<double>(variable_array, pos, time, dt);

        auto dk_rel_ds_L =
            medium
                ->property(MaterialPropertyLib::PropertyType::
                               relative_permeability_nonwetting_phase)
                .template dValue<double>(
                    variable_array,
                    MaterialPropertyLib::Variable::liquid_saturation, pos, time,
                    dt);

        ASSERT_NEAR(k_rel, ref_k_rel_G[idx], 1.0e-10);
        ASSERT_NEAR(dk_rel_ds_L, ref_dk_rel_G_ds_L[idx], 1.0e-10);
    }
}
