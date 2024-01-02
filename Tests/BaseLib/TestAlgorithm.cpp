/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <numeric>
#include <random>
#include <vector>

#include "BaseLib/Algorithm.h"

TEST(BaseLibAlgorithm, testreorderVector)
{
    const std::size_t size = 100;
    std::vector<double> vec(size);
    std::default_random_engine random_engine;
    std::generate(vec.begin(), vec.end(), random_engine);
    std::vector<double> vec0 = vec;

    std::vector<int> order(size);
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), random_engine);

    BaseLib::reorderVector(vec, order);

    for (std::size_t i = 0; i < size; i++)
    {
        EXPECT_EQ(vec0[order[i]], vec[i]);
    }
}

TEST(BaseLibAlgorithm, excludeObjectCopy)
{
    // create and fill vector
    std::size_t const size(100);
    std::vector<std::size_t> v(size);
    std::iota(v.begin(), v.end(), 0);

    std::vector<std::size_t> ex_positions(size / 10);
    // do not copy first 10 elements
    std::iota(ex_positions.begin(), ex_positions.end(), 0);

    std::vector<std::size_t> c1(BaseLib::excludeObjectCopy(v, ex_positions));
    ASSERT_EQ(size - ex_positions.size(), c1.size());

    for (std::size_t i(0); i < c1.size(); i++)
    {
        ASSERT_EQ(c1[i], v[size / 10 + i]);
    }

    // do not copy element 0, 2, 4, 6, 8, 10, 12, 14, 16, 18
    std::transform(ex_positions.begin(), ex_positions.end(),
                   ex_positions.begin(),
                   [](std::size_t const& x) { return x * 2; });

    std::vector<std::size_t> c2(BaseLib::excludeObjectCopy(v, ex_positions));
    ASSERT_EQ(size - ex_positions.size(), c2.size());

    for (std::size_t i(0); i < ex_positions.size(); i++)
    {
        ASSERT_EQ(c2[i], v[2 * i + 1]);
    }
    for (std::size_t i(ex_positions.size()); i < c2.size(); i++)
    {
        ASSERT_EQ(c2[i], v[ex_positions.size() + i]);
    }

    // do not copy the last element
    ex_positions.clear();
    ex_positions.push_back(99);

    std::vector<std::size_t> c3(BaseLib::excludeObjectCopy(v, ex_positions));
    ASSERT_EQ(size - ex_positions.size(), c3.size());

    for (std::size_t i(0); i < c3.size(); i++)
    {
        ASSERT_EQ(c3[i], v[i]);
    }
}

TEST(BaseLibAlgorithm, AnyOf)
{
    namespace BL = BaseLib;

    // empty array
    {
        constexpr std::array<bool, 0> arr = {};
        static_assert(BL::any_of(arr) == false);
    }

    // single element
    {
        constexpr std::array<bool, 1> arr = {false};
        static_assert(BL::any_of(arr) == false);
    }
    {
        constexpr std::array<bool, 1> arr = {true};
        static_assert(BL::any_of(arr) == true);
    }

    // multiple elements
    {
        constexpr std::array<bool, 5> arr = {false, false, false, false, false};
        static_assert(BL::any_of(arr) == false);
    }
    {
        constexpr std::array<bool, 4> arr = {true, true, false, false};
        static_assert(BL::any_of(arr) == true);
    }
    {
        constexpr std::array<bool, 3> arr = {false, true};
        static_assert(BL::any_of(arr) == true);
    }
    {
        constexpr std::array<bool, 5> arr = {false, false, false, false, true};
        static_assert(BL::any_of(arr) == true);
    }
}

TEST(BaseLibAlgorithm, AllOf)
{
    namespace BL = BaseLib;

    // empty array
    {
        constexpr std::array<bool, 0> arr = {};
        static_assert(BL::all_of(arr) == true);
    }

    // single element
    {
        constexpr std::array<bool, 1> arr = {false};
        static_assert(BL::all_of(arr) == false);
    }
    {
        constexpr std::array<bool, 1> arr = {true};
        static_assert(BL::all_of(arr) == true);
    }

    // multiple elements
    {
        constexpr std::array<bool, 4> arr = {true, true, true, true};
        static_assert(BL::all_of(arr) == true);
    }
    {
        constexpr std::array<bool, 5> arr = {false, true, true, true, true};
        static_assert(BL::all_of(arr) == false);
    }
    {
        constexpr std::array<bool, 4> arr = {true, false, true, true};
        static_assert(BL::all_of(arr) == false);
    }
    {
        constexpr std::array<bool, 5> arr = {true, true, true, true, false};
        static_assert(BL::all_of(arr) == false);
    }
    {
        constexpr std::array<bool, 5> arr = {true, false, false, true, false};
        static_assert(BL::all_of(arr) == false);
    }
}

TEST(BaseLibAlgorithm, NoneOf)
{
    namespace BL = BaseLib;

    // empty array
    {
        constexpr std::array<bool, 0> arr = {};
        static_assert(BL::none_of(arr) == true);
    }

    // single element
    {
        constexpr std::array<bool, 1> arr = {false};
        static_assert(BL::none_of(arr) == true);
    }
    {
        constexpr std::array<bool, 1> arr = {true};
        static_assert(BL::none_of(arr) == false);
    }

    // multiple elements
    {
        constexpr std::array<bool, 4> arr = {false, false, false, false};
        static_assert(BL::none_of(arr) == true);
    }
    {
        constexpr std::array<bool, 4> arr = {true, false, false, false};
        static_assert(BL::none_of(arr) == false);
    }
    {
        constexpr std::array<bool, 3> arr = {false, true, false};
        static_assert(BL::none_of(arr) == false);
    }
    {
        constexpr std::array<bool, 5> arr = {false, false, false, false, true};
        static_assert(BL::none_of(arr) == false);
    }
    {
        constexpr std::array<bool, 5> arr = {false, true, false, false, true};
        static_assert(BL::none_of(arr) == false);
    }
}
