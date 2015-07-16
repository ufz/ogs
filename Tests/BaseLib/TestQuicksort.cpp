/**
 * \file
 * \author Dmitrij Naumov
 * \date   Nov. 2012
 * \brief  Testing the specialized version of quicksort.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "gtest/gtest.h"
#include "autocheck/autocheck.hpp"
#include "quicksort.h"

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

namespace ac = autocheck;

struct BaseLibQuicksort : public ::testing::Test
{
    virtual void SetUp()
    {
        cls.trivial([](const std::vector<int>& xs)
                    {
                        return xs.size() < 2;
                    });

        cls.collect([](std::vector<int> const& xs)
                    {
                        return xs.size() < 10 ? "short" : "long";
                    });
    }

    ac::gtest_reporter gtest_reporter;
    ac::classifier<std::vector<int>> cls;
};

// Quicksort result is sorted.
TEST_F(BaseLibQuicksort, SortsAsSTLSort)
{
    cls.classify(
        [](std::vector<int> const& xs)
        {
            return std::is_sorted(xs.begin(), xs.end());
        },
        "sorted");

    auto quicksortSortsAsSTLSort = [](const std::vector<int>& xs) -> bool
    {
        std::vector<std::size_t> perm(xs.size());
        if (!xs.empty())
            BaseLib::quicksort((int*)&(xs[0]), 0, xs.size(), &(perm[0]));
        return std::is_sorted(xs.begin(), xs.end());
    };

    ac::check<std::vector<int>>(quicksortSortsAsSTLSort, 100,
                                ac::make_arbitrary<std::vector<int>>(),
                                gtest_reporter, cls);
}

template <typename T>
struct ordered_unique_list_gen
{
    ac::generator<std::vector<T>> source;
    typedef std::vector<T> result_type;

    std::vector<T> operator()(std::size_t size)
    {
        // Generate double as many elements because many will be discarded by
        // applying unique.
        result_type xs(source(2 * size));
        std::sort(xs.begin(), xs.end());
        auto last = std::unique(xs.begin(), xs.end());
        xs.erase(last, xs.end());
        return std::move(xs);
    }
};

// Permutations of sorted, unique vector remain untouched.
TEST_F(BaseLibQuicksort, ReportCorrectPermutations)
{
    auto gen = ac::make_arbitrary(ordered_unique_list_gen<int>());

    auto quicksortCheckPermutations = [](const std::vector<int>& xs)
    {
        std::vector<std::size_t> perm(xs.size());
        std::iota(perm.begin(), perm.end(), 0);

        BaseLib::quicksort((int*)&(xs[0]), 0, xs.size(), &(perm[0]));

        for (std::size_t i = 0; i < perm.size(); ++i)
            if (perm[i] != i)
            {
                std::cerr << i << " " << perm[i] << "\n";
                return false;
            }
        return true;
    };

    ac::check<std::vector<int>>(quicksortCheckPermutations, 100,
            gen, gtest_reporter, cls);
}

// Permutations of reverse sorted, unique vector is also reversed.
TEST_F(BaseLibQuicksort, ReportCorrectPermutationsReverse)
{
    auto reverse = [](std::vector<int>&& xs, std::size_t) -> std::vector<int>
    {
        std::reverse(xs.begin(), xs.end());
        return xs;
    };

    auto gen = ac::make_arbitrary(ac::map(reverse, ordered_unique_list_gen<int>()));

    auto quicksortCheckPermutations = [](const std::vector<int>& xs)
    {
        std::vector<std::size_t> perm(xs.size());
        std::iota(perm.begin(), perm.end(), 0);

        BaseLib::quicksort((int*)&(xs[0]), 0, xs.size(), &(perm[0]));

        for (std::size_t i = 0; i < perm.size(); ++i)
            if (perm[i] != perm.size() - i - 1) return false;
        return true;
    };

    ac::check<std::vector<int>>(quicksortCheckPermutations, 100,
            gen, gtest_reporter, cls);
}
