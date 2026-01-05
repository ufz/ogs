// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "MathLib/KahanSum.h"

#include <gtest/gtest.h>

#include <range/v3/view/reverse.hpp>

#include "Tests/AutoCheckTools.h"

using namespace MathLib;
namespace ac = autocheck;

struct NonNegativeDoubleListGenerator
{
    using result_type = std::vector<double>;
    result_type operator()(std::size_t const size) const
    {
        result_type r;
        r.reserve(size);
        auto positive_double_gen =
            ac::map(&ac::absoluteValue, ac::generator<double>{});
        std::generate_n(std::back_inserter(r), size,
                        ac::fix(size, positive_double_gen));
        return r;
    }
};

struct MathLibKahanSumPropertyTest : public ::testing::Test
{
    ac::gtest_reporter gtest_reporter;
};

TEST_F(MathLibKahanSumPropertyTest, SortedSequences)
{
    auto property = [](const std::vector<double>& sequence)
    {
        // Sort to get an optimal behaviour for sum computation.
        std::vector<double> sorted{sequence};
        std::sort(sorted.begin(), sorted.end(),
                  [](double const a, double const b)
                  { return std::abs(a) < std::abs(b); });
        KahanSum const a{sorted};

        // Reverse sort yields worst sum behaviour resulting in maximum error.
        KahanSum const b{sorted | ranges::views::reverse};

        if (a() == b())
        {
            return true;
        }
        // Testing the relative error to be within the theoretical limits. Keep
        // in mind, that the sequence is constructed with non-negative numbers
        // and the condition number is 1.
        return std::abs(a() - b()) / ((a() + b()) / 2.) <
               2. * std::numeric_limits<double>::epsilon();
    };

    // The non-negative double list generator is used to keep the condition
    // number, which is defined as sum|a_i|/|sum a_i|, equal to 1.
    autocheck::check<std::vector<double>>(
        property, 1000, ac::make_arbitrary(NonNegativeDoubleListGenerator{}),
        gtest_reporter);
}
