/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   TestLiquidFlowMaterialProperties.cpp
 *
 * Created on August 18, 2016, 4:53 PM
 */

#include <gtest/gtest.h>

#include "TestTools.h"

#include "ProcessLib/LiquidFlow/LiquidFlowMaterialProperties.h"
#include "MaterialLib/Fluid/FluidProperty.h"

#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"

#include "MeshLib/Mesh.h"

using namespace ProcessLib::LiquidFlow;
using namespace MaterialLib::Fluid;
using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;
using Matrix = Eigen::MatrixXd;

TEST(ProcessLibLiquidFlow, checkLiquidFlowMaterialProperties)
{
    const char xml[] =
        "<material_property>"
        "    <fluid>"
        "        <fluid>"
        "            <density>"
        "               <type>LiquidDensity</type>"
        "               <temperature0> 273.15 </temperature0>"
        "               <p0> 1.e+5 </p0> "
        "               <bulk_modulus> 2.15e+9 </bulk_modulus>"
        "               <beta> 2.0e-4 </beta> "
        "               <rho0>999.8</rho0>"
        "            </density>"
        "            <viscosity>"
        "                <type>Vogels</type>"
        "                <liquid_type>Water </liquid_type>"
        "            </viscosity>"
        "        </fluid>"
        "    </fluid>"
        "    <porous_medium>"
        "        <porous_medium>"
        "            <permeability>"
        "                <values>2.e-10 0. 0. 0. 3.e-10 0. 0. 0. 4.0e-10</values>"
        "            </permeability>"
        "            <porosity>"
        "                <type>Constant</type>"
        "                <value> 0.2 </value> "
        "            </porosity>"
        "            <storage>"
        "                <type>Constant</type>"
        "                <value> 1.e-4 </value>"
        "            </storage>"
        "        </porous_medium>"
        "    </porous_medium>"
        "</material_property>";

    auto const ptree = readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("material_property");

    LiquidFlowMaterialProperties lprop(sub_config);

    // Check density
    const ArrayType vars = {273.15 + 60.0, 1.e+6};
    const double T0 = 273.15;
    const double p0 = 1.e+5;
    const double rho0 = 999.8;
    const double K = 2.15e+9;
    const double beta = 2.e-4;
    const double T = vars[0];
    const double p = vars[1];

    const double fac_T = 1. + beta * (T - T0);
    ASSERT_NEAR(rho0 / fac_T / (1. - (p - p0) / K),
                lprop.density_l->getValue(vars), 1.e-10);

    // Test the derivative with respect to temperature.
    ASSERT_NEAR(-beta * rho0 / (fac_T * fac_T) / (1. - (p - p0) / K),
                lprop.density_l->getdValue(vars, PropertyVariableType::T), 1.e-10);

    // Test the derivative with respect to pressure.
    const double fac_p = 1. - (p - p0) / K;
    ASSERT_NEAR(rho0 / (1. + beta * (T - T0)) / (fac_p * fac_p * K),
                lprop.density_l->getdValue(vars, PropertyVariableType::pl), 1.e-10);

    // Check viscosity
    ArrayType vars1;
    vars1[0] = 303.0;
    const auto var_type = MaterialLib::Fluid::PropertyVariableType::T;
    ASSERT_NEAR(0.802657e-3, lprop.viscosity->getValue(vars1), 1.e-5);
    ASSERT_NEAR(-1.87823e-5, lprop.viscosity->getdValue(vars1, var_type),
                1.e-5);

    // Check permeability
    Eigen::MatrixXd& perm = lprop.intrinsic_permeabiliy[0];
    ASSERT_EQ(2.e-10, perm(0, 0));
    ASSERT_EQ(0., perm(0, 1));
    ASSERT_EQ(0., perm(0, 2));
    ASSERT_EQ(0., perm(1, 0));
    ASSERT_EQ(3.e-10, perm(1, 1));
    ASSERT_EQ(0., perm(1, 2));
    ASSERT_EQ(0., perm(2, 0));
    ASSERT_EQ(0., perm(2, 1));
    ASSERT_EQ(4.e-10, perm(2, 2));

    // Check porosity
    const double variable = 0.;
    const double temperature = 0.;
    ASSERT_EQ(0.2, lprop.porosity[0]->getValue(variable, temperature));

    // Check storage
    ASSERT_EQ(1.e-4, lprop.storage[0]->getValue(variable));
}
