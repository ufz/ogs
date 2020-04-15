/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on March 27, 2020, 4:01 PM
 *
 */

#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <limits>
#include <random>

#include "BaseLib/ConfigTree.h"

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/RelativePermeability/CreateRelPermNonWettingPhaseVanGenuchtenMualem.h"
#include "MaterialLib/MPL/Properties/RelativePermeability/RelPermNonWettingPhaseVanGenuchtenMualem.h"

#include "TestMPL.h"
#include "Tests/TestTools.h"

TEST(MaterialPropertyLib, RelPermNonWettingPhaseVanGenuchtenMualem)
{
    const char xml[] =
        "<property>"
        "   <name>relative_permeability</name>"
        "   <type>RelativePermeabilityNonWettingPhaseVanGenuchtenMualem</type>"
        "   <residual_saturation>0.1</residual_saturation>"
        "   <maximum_saturation>1.0</maximum_saturation>"
        "   <exponent>0.5</exponent>"
        "<minimum_relative_permeability>1.e-9</minimum_relative_permeability>"
        "</property>";

    std::unique_ptr<MPL::Property> const perm_ptr = Tests::createTestProperty(
        xml, MPL::createRelPermNonWettingPhaseVanGenuchtenMualem);
    MPL::Property const& perm_model = *perm_ptr;

    std::vector<double> const S = {0.2, 0.33, 0.45, 0.52, 0.6, 0.85};
    std::vector<double> const expected_krel = {
        8.38365641777669e-01, 6.88828521847265e-01, 5.30330085889911e-01,
        4.32869977650940e-01, 3.20750149549792e-01, 2.54616639316144e-02};

    const double perturbation = 1.e-8;
    for (std::size_t i = 0; i < S.size(); i++)
    {
        MPL::VariableArray variable_array;
        ParameterLib::SpatialPosition const pos;
        double const t = std::numeric_limits<double>::quiet_NaN();
        double const dt = std::numeric_limits<double>::quiet_NaN();
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S[i];

        const double k_rel_i =
            perm_model.template value<double>(variable_array, pos, t, dt);

        ASSERT_NEAR(expected_krel[i], k_rel_i, 1.e-9);

        const double dkrel_dS = perm_model.template dValue<double>(
            variable_array, MaterialPropertyLib::Variable::liquid_saturation,
            pos, t, dt);
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S[i] - perturbation;
        const double k_rel_i_0 =
            perm_model.template value<double>(variable_array, pos, t, dt);
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S[i] + perturbation;
        const double k_rel_i_1 =
            perm_model.template value<double>(variable_array, pos, t, dt);

        // Compare the derivative with numerical one.
        ASSERT_NEAR(dkrel_dS, 0.5 * (k_rel_i_1 - k_rel_i_0) / perturbation,
                    1.e-8);
    }
}
