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
#include "MaterialLib/MPL/Properties/AverageMolarMass.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

TEST(MaterialPropertyLib, AverageMolarMass)
{
    std::stringstream m_beg;
    std::stringstream m_end;

    const auto M_Unobtainium = 0.214;
    const auto M_Vibranium = 0.214;

    m_beg << "<medium>\n";
    m_beg << "    <phases>\n";
    m_beg << "        <phase>\n";
    m_beg << "            <type>Gas</type>\n";
    m_beg << "            <components>\n";
    m_beg << "                <component>\n";
    m_beg << "                    <name>Unobtainium</name>\n";
    m_beg << "                    <properties>\n";
    m_beg << "                        <property>\n";
    m_beg << "                            <name>molar_mass</name>\n";
    m_beg << "                            <type>Constant</type>\n";
    m_beg << "                            <value>" << M_Unobtainium
          << "</value>\n";
    m_beg << "                        </property>\n";
    m_beg << "                    </properties>\n";
    m_beg << "                </component>\n";
    m_beg << "                <component>\n";
    m_beg << "                    <name>Vibranium</name>\n";
    m_beg << "                    <properties>\n";
    m_beg << "                        <property>\n";
    m_beg << "                            <name>molar_mass</name>\n";
    m_beg << "                            <type>Constant</type>\n";
    m_beg << "                            <value>" << M_Vibranium
          << "</value>\n";
    m_beg << "                        </property>\n";
    m_beg << "                    </properties>\n";
    m_beg << "                </component>\n";
    m_beg << "            </components>\n";
    m_beg << "            <properties>\n";
    m_beg << "                <property>\n";
    m_beg << "                    <name>mole_fraction</name>\n";
    m_beg << "                    <type>Constant</type>\n";

    m_end << "                </property>\n";
    m_end << "                <property>\n";
    m_end << "                    <name>molar_mass</name>\n";
    m_end << "                    <type>AverageMolarMass</type>\n";
    m_end << "                </property>\n";
    m_end << "            </properties>\n";
    m_end << "        </phase>\n";
    m_end << "    </phases>\n";
    m_end << "</medium>\n";

    std::array<double, 24> ref_molar_fraction_A = {
        0.0000000000000000, 0.0173972718775814, 0.1144192284181090,
        0.1630553923269180, 0.2031396804982220, 0.2518083568978420,
        0.2990064121013330, 0.3275246133944790, 0.3863807424698890,
        0.3944558568846510, 0.4844391099526180, 0.5396876124497050,
        0.5614261169216030, 0.6532570029184130, 0.6881977591185020,
        0.6913221859105890, 0.7309070419102210, 0.8073016060107850,
        0.8487258523827500, 0.9158857855082890, 0.9212465868752240,
        0.9846022661096390, 0.9976737941621020, 1.0000000000000000};

    for (size_t idx = 0; idx < ref_molar_fraction_A.size(); idx++)
    {
        std::stringstream m_composition;

        m_composition << "<value>"
                      << std::setprecision(
                             std::numeric_limits<double>::digits10)
                      << ref_molar_fraction_A[idx] << " "
                      << 1. - ref_molar_fraction_A[idx] << "</value>\n";
        std::stringstream m;
        m << m_beg.str() << m_composition.str() << m_end.str();

        auto const& medium = Tests::createTestMaterial(m.str());
        auto const& gas_phase = medium->phase("Gas");
        MaterialPropertyLib::VariableArray variable_array;
        ParameterLib::SpatialPosition const pos;
        double const time = std::numeric_limits<double>::quiet_NaN();
        double const dt = std::numeric_limits<double>::quiet_NaN();

        auto M =
            gas_phase.property(MaterialPropertyLib::PropertyType::molar_mass)
                .template value<double>(variable_array, pos, time, dt);

        auto const M_ref = M_Unobtainium * ref_molar_fraction_A[idx] +
                           M_Vibranium * (1. - ref_molar_fraction_A[idx]);

        ASSERT_NEAR(M, M_ref, 1.0e-10);
    }
}
