/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on March 27, 2020, 3:01 PM
 *
 */

#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <limits>
#include <random>

#include "BaseLib/ConfigTree.h"

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/CapillaryPressureBrooksCorey.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/CapillaryPressureLiakopoulos.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/CapillaryPressureVanGenuchten.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationBrooksCorey.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationLiakopoulos.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationVanGenuchten.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/CreateCapillaryPressureBrooksCorey.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/CreateCapillaryPressureLiakopoulos.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/CreateCapillaryPressureVanGenuchten.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/CreateCapillaryPressureRegularizedVanGenuchten.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/CreateSaturationBrooksCorey.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/CreateSaturationLiakopoulos.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/CreateSaturationVanGenuchten.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/CapillaryPressureRegularizedVanGenuchten.h"

#include "TestMPL.h"
#include "Tests/TestTools.h"

namespace MPL = MaterialPropertyLib;

class TestCapillaryPressure : public ::testing::Test
{
public:
    void test()
    {
        MPL::Property const& capillary_pressure_property =
            *_capillary_pressure_property;

        MPL::VariableArray variable_array;
        ParameterLib::SpatialPosition const pos;
        double const t = std::numeric_limits<double>::quiet_NaN();
        double const dt = std::numeric_limits<double>::quiet_NaN();

        std::random_device rd;
        std::mt19937 mt(rd());
        const double offset = std::sqrt(std::numeric_limits<double>::epsilon());
        std::uniform_real_distribution<double> distributor(
            _residual_liquid_saturation + offset,
            _maximum_liquid_saturation - offset);

        const int n = 20;
        for (int i = 0; i <= n; ++i)
        {
            double const S = distributor(mt);
            variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
                S;

            const double pc =
                capillary_pressure_property.template value<double>(
                    variable_array, pos, t, dt);

            if (_saturation_property)
            {
                variable_array[static_cast<int>(
                    MPL::Variable::capillary_pressure)] = pc;

                MPL::Property const& saturation_property =
                    *_saturation_property;
                const double computed_saturation =
                    saturation_property.template value<double>(variable_array,
                                                               pos, t, dt);

                ASSERT_LE(std::fabs(S - computed_saturation), 1e-9)
                    << "for saturation " << S
                    << " and re-computed saturation via capillary pressure"
                    << computed_saturation;
            }
            const double dPcdS =
                capillary_pressure_property.template dValue<double>(
                    variable_array,
                    MaterialPropertyLib::Variable::liquid_saturation, pos, t,
                    dt);

            const double d2PcdS2 =
                capillary_pressure_property.template d2Value<double>(
                    variable_array,
                    MaterialPropertyLib::Variable::liquid_saturation,
                    MaterialPropertyLib::Variable::liquid_saturation, pos, t,
                    dt);

            const double S1 =
                std::clamp(S + offset, _residual_liquid_saturation + 2 * offset,
                           _maximum_liquid_saturation - 2 * offset);
            variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
                S1;

            if (std::fabs(S1 - S) < std::numeric_limits<double>::epsilon())
            {
                continue;
            }

            const double pc1 =
                capillary_pressure_property.template value<double>(
                    variable_array, pos, t, dt);
            const double numerical_dPcdS = (pc1 - pc) / (S1 - S);
            ASSERT_LE(std::fabs((dPcdS - numerical_dPcdS) / dPcdS), 1e-4)
                << "for the analytic derivative of dPc/dS " << dPcdS
                << " and numeric derivative of dPc/dS." << numerical_dPcdS;

            const double dPcdS1 =
                capillary_pressure_property.template dValue<double>(
                    variable_array,
                    MaterialPropertyLib::Variable::liquid_saturation, pos, t,
                    dt);

            const double numerical_d2PcdS2 = (dPcdS1 - dPcdS) / (S1 - S);

            ASSERT_LE(std::fabs((d2PcdS2 - numerical_d2PcdS2) / d2PcdS2), 1e-4)
                << "for the analytic second order derivative of d2Pc/dS2 "
                << d2PcdS2
                << " and numeric  second order  derivative of d2Pc/dS2."
                << numerical_d2PcdS2;
        }
    }

protected:
    std::unique_ptr<MaterialPropertyLib::Property> _saturation_property;
    std::unique_ptr<MaterialPropertyLib::Property> _capillary_pressure_property;
    double _residual_liquid_saturation;
    double _maximum_liquid_saturation;
};

std::unique_ptr<MaterialPropertyLib::Property> createTestProperty(
    const char xml[],
    std::function<std::unique_ptr<MaterialPropertyLib::Property>(
        BaseLib::ConfigTree const& config)>
        createProperty)
{
    auto const ptree = readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("property");
    // Parsing the property name:
    auto const property_name =
        sub_config.getConfigParameter<std::string>("name");

    return createProperty(sub_config);
}

TEST_F(TestCapillaryPressure, CapillaryPressureVanGenuchten)
{
    const char xml_S_pc[] =
        "<property>"
        "   <name>saturation</name>"
        "   <type>SaturationVanGenuchten</type>"
        "   <residual_liquid_saturation>0.1</residual_liquid_saturation>"
        "   <residual_gas_saturation>0.05</residual_gas_saturation>"
        "   <exponent>0.79</exponent>"
        "   <p_b>5000.0</p_b>"
        "</property>";
    _saturation_property =
        createTestProperty(xml_S_pc, MPL::createSaturationVanGenuchten);

    const char xml_pc_S[] =
        "<property>"
        "   <name>capillary_pressure</name>"
        "   <type>CapillaryPressureVanGenuchten</type>"
        "   <residual_liquid_saturation>0.1</residual_liquid_saturation>"
        "   <maximum_liquid_saturation>0.95</maximum_liquid_saturation>"
        "   <exponent>0.79</exponent>"
        "   <p_b>5000.0</p_b>"
        "   <maximum_capillary_pressure>1.e+20</maximum_capillary_pressure>"
        "</property>";
    _capillary_pressure_property =
        createTestProperty(xml_pc_S, MPL::createCapillaryPressureVanGenuchten);

    _residual_liquid_saturation = 0.1;
    _maximum_liquid_saturation = 0.95;
    test();
}

TEST_F(TestCapillaryPressure, CapillaryPressureRegularizedVanGenuchten)
{
    const char xml_pc_S[] =
        "<property>"
        "   <name>capillary_pressure</name>"
        "   <type>CapillaryPressureRegularizedVanGenuchten</type>"
        "   <residual_liquid_saturation>0.1</residual_liquid_saturation>"
        "   <maximum_liquid_saturation>0.95</maximum_liquid_saturation>"
        "   <exponent>0.79</exponent>"
        "   <p_b>5000.0</p_b>"
        "</property>";
    _capillary_pressure_property = createTestProperty(
        xml_pc_S, MPL::createCapillaryPressureRegularizedVanGenuchten);

    _residual_liquid_saturation = 0.1;
    _maximum_liquid_saturation = 0.95;
    test();
}

TEST_F(TestCapillaryPressure, CapillaryPressureBrooksCorey)
{
    const char xml_S_pc[] =
        "<property>"
        "   <name>saturation</name>"
        "   <type>SaturationBrooksCorey</type>"
        "   <residual_liquid_saturation>0.12</residual_liquid_saturation>"
        "   <residual_gas_saturation>0.06</residual_gas_saturation>"
        "   <lambda>2.5</lambda>"
        "   <entry_pressure>5678.54</entry_pressure>"
        "</property>";
    _saturation_property =
        createTestProperty(xml_S_pc, MPL::createSaturationBrooksCorey);

    const char xml_pc_S[] =
        "<property>"
        "   <name>capillary_pressure</name>"
        "   <type>CapillaryPressureBrooksCorey</type>"
        "   <residual_liquid_saturation>0.12</residual_liquid_saturation>"
        "   <maximum_liquid_saturation>0.94</maximum_liquid_saturation>"
        "   <lambda>2.5</lambda>"
        "   <entry_pressure>5678.54</entry_pressure>"
        "   <maximum_capillary_pressure>1.e+20</maximum_capillary_pressure>"
        "</property>";
    _capillary_pressure_property =
        createTestProperty(xml_pc_S, MPL::createCapillaryPressureBrooksCorey);

    _residual_liquid_saturation = 0.12;
    _maximum_liquid_saturation = 0.94;
    test();
}

TEST_F(TestCapillaryPressure, CapillaryPressureLiakopoulos)
{
    _saturation_property = std::make_unique<MPL::SaturationLiakopoulos>();
    _capillary_pressure_property =
        std::make_unique<MPL::CapillaryPressureLiakopoulos>();

    _residual_liquid_saturation = 0.2;
    _maximum_liquid_saturation = 1.0;
    test();
}
