/**
 * \file
 * \author Dmitrij Naumov
 * \date   Nov. 2012
 * \brief  Testing the specialized version of quicksort.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <algorithm>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>
#include <autocheck/autocheck.hpp>

#include "BaseLib/quicksort.h"

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

    auto quicksortSortsAsSTLSort = [](std::vector<int>& xs) -> bool
    {
        std::vector<std::size_t> perm(xs.size());
        if (!xs.empty())
            BaseLib::quicksort(xs, 0, xs.size(), perm);
        return std::is_sorted(xs.begin(), xs.end());
    };

    ac::check<std::vector<int>>(quicksortSortsAsSTLSort, 100,
                                ac::make_arbitrary<std::vector<int>>(),
                                gtest_reporter, cls);
}

template <typename T>
struct OrderedUniqueListGen
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

// Permutations of non-empty, sorted, unique vector remain untouched.
TEST_F(BaseLibQuicksort, ReportCorrectPermutations)
{
    auto gen = ac::make_arbitrary(OrderedUniqueListGen<int>());

    auto quicksortCheckPermutations = [](std::vector<int>& xs)
    {
        std::vector<std::size_t> perm(xs.size());
        std::iota(perm.begin(), perm.end(), 0);

        BaseLib::quicksort(xs, 0, xs.size(), perm);

        for (std::size_t i = 0; i < perm.size(); ++i)
            if (perm[i] != i)
            {
                std::cerr << i << " " << perm[i] << "\n";
                return false;
            }
        return true;
    };

    ac::check<std::vector<int>>(quicksortCheckPermutations, 100,
                                gen.discard_if([](std::vector<int> xs)
                                               {
                                                   return xs.empty();
                                               }),
                                gtest_reporter, cls);
}

// Permutations of non-empty, sorted, unique vector remain untouched.
TEST_F(BaseLibQuicksort, ReportCorrectPermutationsWithPointer)
{
    auto gen = ac::make_arbitrary(OrderedUniqueListGen<int>());

    auto quicksortCheckPermutations = [](std::vector<int>& xs)
    {
        std::vector<std::size_t> perm(xs.size());
        std::iota(perm.begin(), perm.end(), 0);

        std::vector<int*> p_xs;
        for (int & x : xs)
            p_xs.push_back(&x);

        BaseLib::quicksort(p_xs, 0, p_xs.size(), perm);

        for (std::size_t i = 0; i < perm.size(); ++i)
            if (perm[i] != i)
            {
                std::cerr << i << " " << perm[i] << "\n";
                return false;
            }
        return true;
    };

    ac::check<std::vector<int>>(quicksortCheckPermutations, 100,
                                gen.discard_if([](std::vector<int> xs)
                                               {
                                                   return xs.empty();
                                               }),
                                gtest_reporter, cls);
}

// Permutations of non-empty, reverse sorted, unique vector is also reversed.
TEST_F(BaseLibQuicksort, ReportCorrectPermutationsReverse)
{
    auto reverse = [](std::vector<int>&& xs, std::size_t) -> std::vector<int>
    {
        std::reverse(xs.begin(), xs.end());
        return xs;
    };

    auto gen =
        ac::make_arbitrary(ac::map(reverse, OrderedUniqueListGen<int>()));

    auto quicksortCheckPermutations = [](std::vector<int>& xs)
    {
        std::vector<std::size_t> perm(xs.size());
        std::iota(perm.begin(), perm.end(), 0);

        BaseLib::quicksort(xs, 0, xs.size(), perm);

        for (std::size_t i = 0; i < perm.size(); ++i)
            if (perm[i] != perm.size() - i - 1) return false;
        return true;
    };

    ac::check<std::vector<int>>(quicksortCheckPermutations, 100,
                                gen.discard_if([](std::vector<int> xs)
                                               {
                                                   return xs.empty();
                                               }),
                                gtest_reporter, cls);
}

// Permutations of non-empty, reverse sorted, unique vector is also reversed.
TEST_F(BaseLibQuicksort, ReportCorrectPermutationsReverseWithPointer)
{
    auto reverse = [](std::vector<int>&& xs, std::size_t) -> std::vector<int>
    {
        std::reverse(xs.begin(), xs.end());
        return xs;
    };

    auto gen =
        ac::make_arbitrary(ac::map(reverse, OrderedUniqueListGen<int>()));

    auto quicksortCheckPermutations = [](std::vector<int>& xs)
    {
        std::vector<std::size_t> perm(xs.size());
        std::iota(perm.begin(), perm.end(), 0);

        std::vector<int*> p_xs;
        for (int & x : xs)
            p_xs.push_back(&x);

        BaseLib::quicksort(p_xs, 0, p_xs.size(), perm);

        for (std::size_t i = 0; i < perm.size(); ++i)
            if (perm[i] != perm.size() - i - 1) return false;
        return true;
    };

    ac::check<std::vector<int>>(quicksortCheckPermutations, 100,
                                gen.discard_if([](std::vector<int> xs)
                                               {
                                                   return xs.empty();
                                               }),
                                gtest_reporter, cls);
}

template <typename T, typename Gen = ac::generator<T>>
struct randomSortedPairGenerator
{
    Gen source;
    using result_type = std::tuple<T, T>;

    result_type operator()(std::size_t)
    {
        auto p = std::make_tuple(source(), source());
        if (std::get<0>(p) > std::get<1>(p))
            std::swap(std::get<0>(p), std::get<1>(p));

        return p;
    }
};

TEST_F(BaseLibQuicksort, SortsRangeAsSTLSort)
{
    auto quicksortSortsRangeAsSTLSort = [](std::vector<int>& xs,
        std::tuple<std::size_t, std::size_t> const& range) -> bool
    {
        if (xs.empty())
            return true;

        std::vector<std::size_t> perm(xs.size());
        std::iota(perm.begin(), perm.end(), 0);
        BaseLib::quicksort(xs, std::get<0>(range), std::get<1>(range), perm);
        return
            std::is_sorted(
                xs.begin() + std::get<0>(range), xs.begin() + std::get<1>(range))
            && std::is_sorted(perm.begin(), perm.begin() + std::get<0>(range))
            && std::is_sorted(perm.begin() + std::get<1>(range), perm.end());
    };

    ac::check<std::vector<int>, std::tuple<std::size_t, std::size_t>>(
        quicksortSortsRangeAsSTLSort,
        100,
        ac::make_arbitrary(ac::generator<std::vector<int>>(),
                           randomSortedPairGenerator<std::size_t>())
            .discard_if([](std::vector<int> const& xs,
                           std::tuple<std::size_t, std::size_t> const& range)
                        {
                            return std::get<1>(range) >= xs.size();
                        }),
        gtest_reporter);
}
