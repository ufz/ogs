/**
 *  \copyright
 *   Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   TestMathTools.cpp
 *
 * Created on December 1, 2016, 3:06 PM
 */

#include <gtest/gtest.h>

#include "MathLib/MathTools.h"

using namespace MathLib;

TEST(MathLibMathTools, check_integer_pow)
{
    const double base = 5.0;
    std::size_t n = 10;
    double expected_result = base;
    for (std::size_t i = 1; i < n; i++)
    {
        ASSERT_EQ(expected_result, integer_pow(base, i));
        expected_result *= base;
    }
}
