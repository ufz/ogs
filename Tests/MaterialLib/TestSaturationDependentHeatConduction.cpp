/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on November 10, 2020, 10:34 AM
 */

#include <gtest/gtest.h>

#include <Eigen/Eigen>
#include <boost/math/constants/constants.hpp>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/CreateSaturationDependentHeatConduction.h"
#include "MaterialLib/MPL/Properties/SaturationDependentHeatConduction.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "ParameterLib/ConstantParameter.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

TEST(MaterialPropertyLib, SaturationDependentHeatConduction)
{
    double const k_dry = 0.2;
    double const k_wet = 1.5;

    auto const k_model_eff = MPL::SaturationDependentHeatConduction(
            "thermal_conductivity", k_dry, k_wet);


    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    MPL::VariableArray vars;

    /// Just test for an arbitrary set of values
    {
        double const sat = 0.345;

        vars[static_cast<int>(MPL::Variable::liquid_saturation)] = sat;

        auto const k = std::get<double>(k_model_eff.value(vars, pos, t, dt));

        double const k_expected = 0.6485;

        ASSERT_LE(std::fabs(k_expected - k) / k_expected, 1e-10)
            << "for expected effective thermal conductivity " << k_expected
            << " and for computed thermal conductivity " << k;
    }
}
