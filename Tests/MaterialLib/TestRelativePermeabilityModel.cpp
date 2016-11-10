/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   TestRelativePermeabilityModel.cpp
 *
 * Created on November 2, 2016, 3:09 PM
 */

#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "BaseLib/ConfigTree.h"

#include "TestTools.h"

#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/RelativePermeability.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/CreateRelativePermeabilityModel.h"

using namespace MaterialLib;
using namespace MaterialLib::PorousMedium;

std::unique_ptr<RelativePermeability> createRelativePermeabilityModel(
    const char xml[])
{
    auto const ptree = readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("relative_permeability");
    sub_config.ignoreConfigAttribute("id");
    return MaterialLib::PorousMedium::createRelativePermeabilityModel(
        sub_config);
}

TEST(MaterialPorousMedium, checkWettingPhaseVanGenuchten)
{
    const char xml[] =
        // Should be <relative_permeability id="0">
        "<relative_permeability id=\"0\">"
        "   <type>WettingPhaseVanGenuchten</type>"
        "   <sr>  0.2772 </sr> "
        "   <smax> 1. </smax> "
        "   <m> 0.5 </m> "
        "   <krel_min> 1.e-9 </krel_min> "
        "</relative_permeability>";
    auto const perm_model = createRelativePermeabilityModel(xml);

    std::vector<double> S = {0.2, 0.33, 0.45, 0.52, 0.6, 0.85};
    std::vector<double> krel = {1.e-9,
                                1.9291770102827742e-006,
                                0.00041114113180510661,
                                0.0019569840357319015,
                                0.0074049488610241927,
                                0.13546615958442240};

    for (std::size_t i = 0; i < S.size(); i++)
    {
        ASSERT_NEAR(krel[i], perm_model->getValue(S[i]), 1.e-9);
    }
}

TEST(MaterialPorousMedium, checkNonWettingPhaseVanGenuchten)
{
    const char xml[] =
        "<relative_permeability id=\"0\">"
        "   <type>NonWettingPhaseVanGenuchten</type>"
        "   <sr>  0.1 </sr> "
        "   <smax> 1. </smax> "
        "   <m> 0.5 </m> "
        "   <krel_min> 1.e-9 </krel_min> "
        "</relative_permeability>";
    auto const perm_model = createRelativePermeabilityModel(xml);

    std::vector<double> S = {0.2, 0.33, 0.45, 0.52, 0.6, 0.85};
    std::vector<double> krel = {0.87422700239237161, 0.74331414436457388,
                                0.59527539448807487, 0.49976666464188485,
                                0.38520070797257489, 0.041219134248319585};

    for (std::size_t i = 0; i < S.size(); i++)
    {
        ASSERT_NEAR(krel[i], perm_model->getValue(S[i]), 1.e-9);
    }
}

TEST(MaterialPorousMedium, checkWettingPhaseBrookCoreyOilGas)
{
    const char xml[] =
        "<relative_permeability id=\"0\">"
        "   <type>WettingPhaseBrookCoreyOilGas</type>"
        "   <sr>  0.2 </sr> "
        "   <smax> 0.8 </smax> "
        "   <m> 2 </m> "
        "   <krel_min> 1.e-9 </krel_min> "
        "</relative_permeability>";
    auto const perm_model = createRelativePermeabilityModel(xml);

    std::vector<double> S = {0.2, 0.33, 0.45, 0.52, 0.6, 0.85};
    std::vector<double> krel = {0.,
                                0.0022037808641975302,
                                0.030140817901234556,
                                0.080908641975308573,
                                0.19753086419753069,
                                1.0};

    for (std::size_t i = 0; i < S.size(); i++)
    {
        ASSERT_NEAR(krel[i], perm_model->getValue(S[i]), 1.e-9);
    }
}

TEST(MaterialPorousMedium, checkNonWettingPhaseBrookCoreyOilGas)
{
    const char xml[] =
        "<relative_permeability id=\"0\">"
        "   <type>NonWettingPhaseBrookCoreyOilGas</type>"
        "   <sr>  0.2 </sr> "
        "   <smax> 0.8 </smax> "
        "   <m> 2 </m> "
        "   <krel_min> 1.e-9 </krel_min> "
        "</relative_permeability>";
    auto const perm_model = createRelativePermeabilityModel(xml);

    std::vector<double> S = {0.2, 0.33, 0.45, 0.52, 0.6, 0.85};
    std::vector<double> krel = {1.0,
                                0.58480547839506147,
                                0.28120177469135793,
                                0.15583209876543211,
                                0.061728395061728412,
                                .0};

    for (std::size_t i = 0; i < S.size(); i++)
    {
        ASSERT_NEAR(krel[i], perm_model->getValue(S[i]), 1.e-9);
    }
}

TEST(MaterialPorousMedium, checkReletivePermeabilityCurve)
{
    const char xml[] =
        "<relative_permeability id=\"0\">"
        "   <type>Curve</type>"
        "       <curve>"
        "           <point><data> 0. 0.9   </data></point>"
        "           <point><data> 0.4 0.5  </data></point>"
        "           <point><data> 0.9 0.01 </data></point>"
        "       </curve>"
        "</relative_permeability>";
    auto const perm_model = createRelativePermeabilityModel(xml);

    std::vector<double> S = {0.2, 0.33, 0.45, 0.52, 0.6, 0.85};
    std::vector<double> krel = {0.7, 0.57, 0.451, 0.3824, 0.304, 0.059};

    for (std::size_t i = 0; i < S.size(); i++)
    {
        ASSERT_NEAR(krel[i], perm_model->getValue(S[i]), 1.e-9);
    }
}
