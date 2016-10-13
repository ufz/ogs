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
#include <memory>

#include "TestTools.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"

#include "ProcessLib/Parameter/SpatialPosition.h"
#include "ProcessLib/LiquidFlow/LiquidFlowMaterialProperties.h"

#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"

using namespace ProcessLib::LiquidFlow;
using namespace MaterialLib::Fluid;
using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;
using Matrix = Eigen::MatrixXd;

TEST(ProcessLibLiquidFlow, checkLiquidFlowMaterialProperties)
{
    const char xml[] =
        "<material_property>"
        "    <fluid>"
        "        <density>"
        "           <type>LiquidDensity</type>"
        "           <temperature0> 273.15 </temperature0>"
        "           <p0> 1.e+5 </p0> "
        "           <bulk_modulus> 2.15e+9 </bulk_modulus>"
        "           <beta> 2.0e-4 </beta> "
        "           <rho0>999.8</rho0>"
        "        </density>"
        "        <viscosity>"
        "            <type>Vogels</type>"
        "            <liquid_type>Water </liquid_type>"
        "        </viscosity>"
        "    </fluid>"
        "    <porous_medium>"
        "        <porous_medium  id=\"0\">"
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

    MeshLib::Properties dummy_property;
    auto const& dummy_property_vector =
                dummy_property.createNewPropertyVector<int>(
                               "MaterialIDs", MeshLib::MeshItemType::Cell, 1);

    LiquidFlowMaterialProperties lprop(sub_config, *dummy_property_vector);

    ProcessLib::SpatialPosition pos;
    pos.setElementID(0);

    // Check permeability
    const Eigen::MatrixXd& perm = lprop.getPermeability(0., pos, 1);
    ASSERT_EQ(2.e-10, perm(0, 0));
    ASSERT_EQ(0., perm(0, 1));
    ASSERT_EQ(0., perm(0, 2));
    ASSERT_EQ(0., perm(1, 0));
    ASSERT_EQ(3.e-10, perm(1, 1));
    ASSERT_EQ(0., perm(1, 2));
    ASSERT_EQ(0., perm(2, 0));
    ASSERT_EQ(0., perm(2, 1));
    ASSERT_EQ(4.e-10, perm(2, 2));

    const double T = 273.15 + 60.0;
    const double p = 1.e+6;
    const double mass_coef = lprop.getMassCoefficient(0., pos, p, T, 0., 0.);
    ASSERT_NEAR(0.000100000093, mass_coef, 1.e-10);
}
