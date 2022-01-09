/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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
        "   <residual_liquid_saturation>0.1</residual_liquid_saturation>"
        "   <residual_gas_saturation>0.05</residual_gas_saturation>"
        "   <exponent>0.5</exponent>"
        "   <min_relative_permeability>1.e-9</min_relative_permeability>"
        "</property>";

    std::unique_ptr<MPL::Property> const perm_ptr = Tests::createTestProperty(
        xml, MPL::createRelPermNonWettingPhaseVanGenuchtenMualem);
    MPL::Property const& perm_model = *perm_ptr;

    // Saturation values including S=0 < S_r, S=1>S_max, which are used to check
    // the saturation restriction, i,e S in [Sr, S_max].
    std::vector<double> const S_L = {
        -1.0e-4, 0.1 - 1.e-9,  0.1,  0.1 + 1.e-9,  0.2, 0.33, 0.45, 0.52, 0.6,
        0.85,    0.95 - 1.e-9, 0.95, 0.95 + 1.e-9, 1.0, 1.1};
    std::vector<double> const expected_krel = {1.0,
                                               1.0,
                                               1.0,
                                               1.0,
                                               9.263352402730153e-01,
                                               7.915237953228073e-01,
                                               6.369259422953941e-01,
                                               5.375997895545541e-01,
                                               4.196512496776206e-01,
                                               7.595785085896588e-02,
                                               1.0e-9,
                                               1.0e-9,
                                               1.0e-9,
                                               1.0e-9,
                                               1.0e-9};

    const double S_L_r = 0.1;     // 1.0 - S_n_max;
    const double S_L_max = 0.95;  // 1.0 - S_n_r;
    for (std::size_t i = 0; i < S_L.size(); i++)
    {
        MPL::VariableArray variable_array;
        ParameterLib::SpatialPosition const pos;
        double const t = std::numeric_limits<double>::quiet_NaN();
        double const dt = std::numeric_limits<double>::quiet_NaN();
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L[i];

        const double k_rel_i =
            perm_model.template value<double>(variable_array, pos, t, dt);

        ASSERT_NEAR(expected_krel[i], k_rel_i, 1.e-9);

        const double dkrel_dS = perm_model.template dValue<double>(
            variable_array, MaterialPropertyLib::Variable::liquid_saturation,
            pos, t, dt);

        double S_L_a = S_L[i];
        double S_L_b = S_L[i];
        double factor = 1.0;
        double perturbation = 1.e-8;

        if (std::fabs(S_L[i] - S_L_r) < std::numeric_limits<double>::epsilon())
        {
            S_L_b += perturbation;
        }
        else if (std::fabs(S_L[i] - S_L_max) <
                 std::numeric_limits<double>::epsilon())
        {
            S_L_a -= perturbation;
        }
        else
        {
            const double S_L = S_L_a;
            S_L_a -= perturbation;
            S_L_b += perturbation;
            if (S_L < S_L_r && S_L_b > S_L_r)
            {
                perturbation = 0.9 * (S_L_r - S_L);
                S_L_a = S_L - perturbation;
                S_L_b = S_L + perturbation;
            }
            else if (S_L > S_L_r && S_L_a < S_L_r)
            {
                perturbation = 0.9 * (S_L - S_L_r);
                S_L_a = S_L - perturbation;
                S_L_b = S_L + perturbation;
            }
            else if (S_L > S_L_max && S_L_a < S_L_max)
            {
                perturbation = 0.9 * (S_L - S_L_max);
                S_L_a = S_L - perturbation;
                S_L_b = S_L + perturbation;
            }
            else if (S_L < S_L_max && S_L_b > S_L_max)
            {
                perturbation = 0.9 * (S_L_max - S_L);
                S_L_a = S_L - perturbation;
                S_L_b = S_L + perturbation;
            }

            factor = 0.5;
        }

        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L_a;
        const double k_rel_i_0 =
            perm_model.template value<double>(variable_array, pos, t, dt);
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L_b;
        const double k_rel_i_1 =
            perm_model.template value<double>(variable_array, pos, t, dt);

        // Compare the derivative with numerical one.
        ASSERT_NEAR(dkrel_dS, factor * (k_rel_i_1 - k_rel_i_0) / perturbation,
                    6.0e-8);
    }

    // Test the calculation of the liquid saturation for the minimum relative
    // permeability.
    MPL::RelPermNonWettingPhaseVanGenuchtenMualem k_non_wetting(
        "dummy", 0.1, 0.4, 0.3288590604, 1.e-9);
    ASSERT_NEAR(0.59999999552634464,
                k_non_wetting.computeSaturationForMinimumRelativePermeability(),
                1.0e-16);
}
