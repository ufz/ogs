/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <gtest/gtest.h>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationVolStrain.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

namespace MPL = MaterialPropertyLib;

TEST(MaterialPropertyLib, SaturationVolStrain)
{
    double const residual_liquid_saturation = 0.1;
    double const residual_gas_saturation = 0.05;
    double const exponent = 0.79;
    double const p_b = 5000;
    double const b11 = 300;
    double const b22 = 100;
    double const b33 = 50;

    MPL::Property const& pressure_saturation =
        MPL::SaturationVolStrain{"saturation", residual_liquid_saturation,
                                    residual_gas_saturation, exponent, p_b, b11, b22, b33};

    MPL::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    double const vol_s = 0.0001;
    double const p_L = 10000;

    variable_array.capillary_pressure = p_L;
    variable_array.volumetric_strain = vol_s;

    double const S = pressure_saturation.template value<double>(
        variable_array, pos, t, dt);

    double const S_expected = 0.158397797;

    ASSERT_NEAR(S, S_expected, 1e-5);

}
