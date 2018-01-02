/**
 * \date   2015-03-18
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <algorithm>
#include <cstddef>
#include <functional>
#include <numeric>
#include <vector>

#include "gtest/gtest.h"

#include "BaseLib/excludeObjectCopy.h"

TEST(BaseLib, excludeObjectCopy)
{
    // create and fill vector
    std::size_t const size(100);
    std::vector<std::size_t> v(size);
    std::iota(v.begin(), v.end(), 0);

    std::vector<std::size_t> ex_positions(size/10);
    // do not copy first 10 elements
    std::iota(ex_positions.begin(), ex_positions.end(), 0);

    std::vector<std::size_t> c1(BaseLib::excludeObjectCopy(v,ex_positions));
    ASSERT_EQ(size-ex_positions.size(), c1.size());

    for (std::size_t i(0); i<c1.size(); i++)
        ASSERT_EQ(c1[i], v[size/10+i]);

    // do not copy element 0, 2, 4, 6, 8, 10, 12, 14, 16, 18
    std::transform(ex_positions.begin(), ex_positions.end(),
        ex_positions.begin(), std::bind1st(std::multiplies<std::size_t>(),2));

    std::vector<std::size_t> c2(BaseLib::excludeObjectCopy(v,ex_positions));
    ASSERT_EQ(size-ex_positions.size(), c2.size());

    for (std::size_t i(0); i<ex_positions.size(); i++)
        ASSERT_EQ(c2[i], v[2*i+1]);
    for (std::size_t i(ex_positions.size()); i<c2.size(); i++)
        ASSERT_EQ(c2[i], v[ex_positions.size()+i]);

    // do not copy the last element
    ex_positions.clear();
    ex_positions.push_back(99);

    std::vector<std::size_t> c3(BaseLib::excludeObjectCopy(v,ex_positions));
    ASSERT_EQ(size-ex_positions.size(), c3.size());

    for (std::size_t i(0); i<c3.size(); i++)
        ASSERT_EQ(c3[i], v[i]);
}
