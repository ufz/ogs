/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *   \file TestreorderVector.cpp
 *
 *   Created on August 14, 2016, 10:01 AM
 *
 */

#include <vector>
#include <algorithm>
#include <numeric>

#include <gtest/gtest.h>

#include "BaseLib/reorderVector.h"

TEST(BaseLib_reorderVector, testreorderVector)
{
    const std::size_t size = 100;
    std::vector<double> vec(size);
    std::generate(vec.begin(), vec.end(), std::rand);
    std::vector<double> vec0 = vec;

    std::vector<int> order(size);
    std::iota(order.begin(), order.end(), 0);
    std::random_shuffle(order.begin(), order.end());

    BaseLib::reorderVector(vec, order);

    for (std::size_t i= 0; i<size; i++)
    {
        EXPECT_EQ(vec0[order[i]], vec[i]);
    }
}
