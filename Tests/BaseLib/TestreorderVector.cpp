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

#include <gtest/gtest.h>

#include "BaseLib/reorderVector.h"

TEST(BaseLib_reorderVector, testreorderVector)
{
    std::vector<double> vec {2016.0, 1996.0, 2006.0, 1986.0};
    std::vector<int> order {3, 1, 2, 0};
    BaseLib::reorderVector(vec, order);

    EXPECT_EQ(1986.0, vec[0]);
    EXPECT_EQ(1996.0, vec[1]);
    EXPECT_EQ(2006.0, vec[2]);
    EXPECT_EQ(2016.0, vec[3]);
}
