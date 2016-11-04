/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   TestCapillaryPressureSaturationModel.cpp
 *
 * Created on November 1, 2016, 11:06 AM
 */

#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "BaseLib/ConfigTree.h"

#include "TestTools.h"

#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CapillaryPressureSaturation.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CreateCapillaryPressureModel.h"

using namespace MaterialLib;
using namespace MaterialLib::PorousMedium;

std::unique_ptr<CapillaryPressureSaturation> createCapillaryPressureModel(
    const char xml[])
{
    auto const ptree = readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("capillary_pressure");
    sub_config.ignoreConfigAttribute("id");
    return MaterialLib::PorousMedium::createCapillaryPressureModel(sub_config);
}

TEST(MaterialPorousMedium, checkBrookCoreyCapillaryPressure)
{
    const char xml[] =
        "<capillary_pressure>"
        "   <type>BrookCorey</type>"
        "   <pd> 19600.0 </pd> "
        "   <sr> 0.1 </sr> "
        "   <smax> 1. </smax> "
        "   <m> 2 </m> "
        "   <pc_max> 1.e11 </pc_max> "
        "</capillary_pressure>";
    auto const pc_model = createCapillaryPressureModel(xml);

    std::vector<double> S = {0.11, 0.2, 0.3, 0.44, 0.52, 0.6, 0.85};
    // The following expected data are calculated by using  OGS5.
    std::vector<double> pc = {185941.926417901, 58800,
                              41577.878733769,  31888.7772993424,
                              28691.4621446869, 26296.1594153975,
                              21470.7242542025};
    std::vector<double> dpc_dS = {-9297096.32089506, -294000,
                                  -103944.696834422, -46895.2607343271,
                                  -34156.5025531986, -26296.1594153975,
                                  -14313.8161694683};

    const double tol_pc = 1.e-6;
    for (std::size_t i = 0; i < S.size(); i++)
    {
        ASSERT_NEAR(S[i], pc_model->getSaturation(pc[i]), 1e-14);
        ASSERT_NEAR(pc[i], pc_model->getCapillaryPressure(S[i]), tol_pc);
        ASSERT_NEAR(dpc_dS[i], pc_model->getdPcdS(S[i]), tol_pc);
    }
}

// In the following test, the expected results are calculated by using OGS5.
TEST(MaterialPorousMedium, checkVanGenuchtenCapillaryPressure)
{
    // rho=1000, alpha = 0.37, pd=rho*g/alpha
    const char xml[] =
        "<capillary_pressure>"
        "   <type>vanGenuchten</type>"
        "   <pd> 26513.513513513513 </pd> "
        "   <sr> 0. </sr> "
        "   <smax> 1. </smax> "
        "   <m> 0.7 </m> "
        "   <pc_max> 1.e5 </pc_max> "
        "</capillary_pressure>";
    auto const pc_model = createCapillaryPressureModel(xml);

    std::vector<double> S = {0.11, 0.2, 0.3, 0.44, 0.52, 0.6, 0.85};
    // The following expected data are calculated by using  OGS5.
    std::vector<double> pc = {
        67392.5090996789, 51197.5842770154, 41864.8480636163, 33730.6152475992,
        30210.4060771713, 27091.7425032625, 17726.625235802};
    std::vector<double> dpc_dS = {-274283.722262232, -121944.994573743,
                                  -72852.9253517949, -47580.0457562877,
                                  -41013.0131889419, -37359.7105247455,
                                  -43138.4488851645};

    const double tol_pc = 1.e-6;
    for (std::size_t i = 0; i < S.size(); i++)
    {
        ASSERT_NEAR(S[i], pc_model->getSaturation(pc[i]), 1e-14);
        ASSERT_NEAR(pc[i], pc_model->getCapillaryPressure(S[i]), tol_pc);
        ASSERT_NEAR(dpc_dS[i], pc_model->getdPcdS(S[i]), tol_pc);
    }
}

TEST(MaterialPorousMedium, checkCapillaryPressureCurve)
{
    const char xml[] =
        "<capillary_pressure id=\"0\">"
        "   <type>Curve</type>"
        "       <curve>"
        "           <coords> 0.     0.5   0.9  </coords>"
        "           <values> 1.e+6  1.e+4 0.0  </values>"
        "       </curve>"
        "</relative_permeability>";
    auto const pc_model = createCapillaryPressureModel(xml);

    std::vector<double> S = {0.0, 0.2, 0.3, 0.52, 0.6, 1.0};
    // The following expected data are calculated by using  OGS5.
    std::vector<double> pc = {1.e+6, 604000.0, 406000.0, 9500.0, 7500.0, 0.};
    std::vector<double> dpc_dS = {-1980000.0, -1980000.0, -1980000.0,
                                  -25000.0,   -25000.0,   -25000.0};

    for (std::size_t i = 0; i < S.size(); i++)
    {
        ASSERT_NEAR(pc[i], pc_model->getCapillaryPressure(S[i]), 1.e-5);
        ASSERT_NEAR(dpc_dS[i], pc_model->getdPcdS(S[i]), 1.e-5);

        if (i == S.size() - 1)
            // Minimum Pc to maximum saturation
            ASSERT_NEAR(0.9, pc_model->getSaturation(pc[i]), 1.e-5);
        else
            ASSERT_NEAR(S[i], pc_model->getSaturation(pc[i]), 1.e-5);
    }
}