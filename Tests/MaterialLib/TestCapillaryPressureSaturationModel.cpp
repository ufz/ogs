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
    std::vector<double> pc = {185941.926417901, 58800,
                              41577.878733769,  31888.7772993424,
                              28691.4621446869, 26296.1594153975,
                              21470.7242542025};
    std::vector<double> dpc_dS = {-9297096.32089506, -294000,
                                  -103944.696834422, -46895.2607343271,
                                  -34156.5025531986, -26296.1594153975,
                                  -14313.8161694683};

    for (std::size_t i = 0; i < S.size(); i++)
    {
        ASSERT_NEAR(pc[i], pc_model->getCapillaryPressure(S[i]), 1.e-5);
        ASSERT_NEAR(S[i], pc_model->getSturation(pc[i]), 1.e-5);
        ASSERT_NEAR(dpc_dS[i], pc_model->getdPcdS(S[i]), 1.e-5);
    }
}

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
    std::vector<double> pc = {
        67392.5090996789, 51197.5842770154, 41864.8480636163, 33730.6152475992,
        30210.4060771713, 27091.7425032625, 17726.625235802};
    std::vector<double> dpc_dS = {-274283.722262232, -121944.994573743,
                                  -72852.9253517949, -47580.0457562877,
                                  -41013.0131889419, -37359.7105247455,
                                  -43138.4488851645};

    for (std::size_t i = 0; i < S.size(); i++)
    {
        ASSERT_NEAR(pc[i], pc_model->getCapillaryPressure(S[i]), 1.e-5);
        ASSERT_NEAR(S[i], pc_model->getSturation(pc[i]), 1.e-5);
        ASSERT_NEAR(dpc_dS[i], pc_model->getdPcdS(S[i]), 1.e-5);
    }
}
