/*!
   \file
   \brief Test class for cubic law permeability model.

   \copyright
    Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#include <gtest/gtest.h>

#include <memory>

#include "MaterialLib/MPL/Properties/CubicLawPermeability.h"
#include "ParameterLib/ConstantParameter.h"
#include "TestMPL.h"

TEST(MaterialPropertyLib, checkCubicLawPermeability)
{
    ParameterLib::ConstantParameter<double> const b("b", 1.e-4);

    MPL::Property const& k_model =
        MPL::CubicLawPermeability{"CubicLawPermeability", b};

    double const expected_k = 8.3333e-10;

    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    MPL::VariableArray vars;

    double const k = k_model.template value<double>(vars, pos, t, dt);

    ASSERT_NEAR(expected_k, k, 3.4e-15);
}
