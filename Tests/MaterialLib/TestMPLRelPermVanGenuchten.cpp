/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <gtest/gtest.h>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/RelativePermeability/RelPermVanGenuchten.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

namespace MPL = MaterialPropertyLib;

TEST(MaterialPropertyLib, RelPermVanGenuchten)
{
    double const residual_liquid_saturation = 0.1;
    double const residual_gas_saturation = 0.05;
    double const exponent = 0.79;
    double const min_relative_permeability_liquid = 1e-12;

    MPL::Property const& permeability = MPL::RelPermVanGenuchten{
        "relative_permeability", residual_liquid_saturation,
        residual_gas_saturation, min_relative_permeability_liquid, exponent};

    MPL::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    double const S_0 = -0.1;
    double const S_max = 1.1;
    int const n_steps = 1000;
    for (int i = 0; i <= n_steps; ++i)
    {
        double const S_L = S_0 + i * (S_max - S_0) / n_steps;
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L;

        double const k_rel =
            permeability.template value<double>(variable_array, pos, t, dt);
        double const dk_rel = permeability.template dValue<double>(
            variable_array, MPL::Variable::liquid_saturation, pos, t, dt);

        double const eps = 1e-8;
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L - eps;
        double const k_rel_minus =
            permeability.template value<double>(variable_array, pos, t, dt);
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L + eps;
        double const k_rel_plus =
            permeability.template value<double>(variable_array, pos, t, dt);

        double const Dk_rel = (k_rel_plus - k_rel_minus) / 2 / eps;

        if (S_L < 1 - residual_gas_saturation)
        {
            ASSERT_LE(std::abs(dk_rel - Dk_rel), 1e-7)
                << "for saturation " << S_L << " and relative permeability "
                << k_rel;
        }
        else
        {
            ASSERT_EQ(dk_rel, 0.) << "for saturation " << S_L
                                  << " and relative permeability " << k_rel;
        }
    }
}
