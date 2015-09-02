/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <gtest/gtest.h>

#include <array>

#include "NumLib/Function/Interpolation.h"

TEST(NumLibFunctionInterpolationTest, TwoVariablesTwoNodes)
{
    double v1 = 0.0, v2 = 0.0;
    std::array<double*, 2> interpolated_values = { &v1, &v2 };

    const std::array<double, 4> nodal_values = {
        0.0, 1.0, // v1
        -1.0, 1.0  // v2
    };

    const std::array<double, 2> shape_matrix = { 0.25, 0.75 };

    NumLib::shapeFunctionInterpolate(nodal_values, shape_matrix,
                                     interpolated_values);

    ASSERT_EQ(0.75, v1);
    ASSERT_EQ(0.5,  v2);
}

