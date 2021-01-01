/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <algorithm>
#include <array>
#include <cstring>
#include <functional>
#include <limits>
#include <random>
#include <string>

#include <gtest/gtest.h>

#include "BaseLib/FileTools.h"

// https://stackoverflow.com/questions/440133/how-do-i-create-a-random-alpha-numeric-string-in-c
template <typename T = std::mt19937>
auto random_generator() -> T
{
    auto constexpr seed_bits = sizeof(typename T::result_type) * T::state_size;
    auto constexpr seed_len =
        seed_bits / std::numeric_limits<std::seed_seq::result_type>::digits;
    auto seed = std::array<std::seed_seq::result_type, seed_len>{};
    auto dev = std::random_device{};
    std::generate_n(begin(seed), seed_len, std::ref(dev));
    auto seed_seq = std::seed_seq(begin(seed), end(seed));
    return T{seed_seq};
}

auto generate_random_alphanumeric_string(std::size_t len) -> std::string
{
    static constexpr auto chars =
        "123456788"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
    thread_local auto rng = random_generator<>();

    auto dist = std::uniform_int_distribution{{}, std::strlen(chars)};
    auto result = std::string(len, '\0');
    std::generate_n(begin(result), len, [&]() { return chars[dist(rng)]; });
    return result;
}

TEST(BaseLib, getParenthesizedString)
{
    {  // {} return empty string
        std::string const test_string = "{}";
        ASSERT_TRUE(std::get<0>(
            BaseLib::getParenthesizedString(test_string, '{', '}', 0)).empty());
    }

    {  // a{}b return empty string
        std::string const pre = generate_random_alphanumeric_string(20);
        std::string const post = generate_random_alphanumeric_string(20);
        std::string const test_string = pre + "{}" + post;
        ASSERT_TRUE(std::get<0>(
            BaseLib::getParenthesizedString(test_string, '{', '}', 0)).empty());
    }

    // }a{
    for (int i = 0; i < 100; ++i)
    {
        std::string const random = generate_random_alphanumeric_string(20);
        std::string const test_string = "}" + random + "{";
        ASSERT_TRUE(std::get<0>(
            BaseLib::getParenthesizedString(test_string, '{', '}', 0)).empty());
    }

    // {a}
    for (int i = 0; i < 100; ++i)
    {
        std::string const expected = generate_random_alphanumeric_string(20);
        std::string const test_string = "{" + expected + "}";
        ASSERT_EQ(expected,
                  std::get<0>(BaseLib::getParenthesizedString(test_string, '{',
                                                              '}', 0)));
        ASSERT_TRUE(std::get<0>(
            BaseLib::getParenthesizedString(test_string, '{', '}', 1)).empty());
    }

    // a{b}
    for (int i = 0; i < 100; ++i)
    {
        std::string const pre = generate_random_alphanumeric_string(20);
        std::string const expected = generate_random_alphanumeric_string(20);
        std::string const test_string = pre + "{" + expected + "}";
        ASSERT_EQ(expected,
                  std::get<0>(BaseLib::getParenthesizedString(test_string, '{',
                                                              '}', 0)));
        ASSERT_TRUE(std::get<0>(BaseLib::getParenthesizedString(
                                    test_string, '{', '}', pre.length() + 2))
                        .empty());
    }

    // a{b}c"
    for (int i = 0; i < 100; ++i)
    {
        std::string const pre = generate_random_alphanumeric_string(20);
        std::string const post = generate_random_alphanumeric_string(20);
        std::string const expected = generate_random_alphanumeric_string(20);
        std::string const test_string = pre + "{" + expected + "}" + post;
        ASSERT_EQ(expected,
                  std::get<0>(BaseLib::getParenthesizedString(test_string, '{',
                                                              '}', 0)));
        ASSERT_TRUE(std::get<0>(
            BaseLib::getParenthesizedString(
                test_string, '{', '}', pre.length() + expected.length() + 2))
                .empty());
    }

    // a{b{c}"
    for (int i = 0; i < 100; ++i)
    {
        std::string const pre = generate_random_alphanumeric_string(20);
        std::string const post = generate_random_alphanumeric_string(20);
        std::string const expected =
            generate_random_alphanumeric_string(20) + "{" + post;
        std::string const test_string = pre + "{" + expected + "}";
        ASSERT_EQ(expected,
                  std::get<0>(BaseLib::getParenthesizedString(test_string, '{',
                                                              '}', 0)));
        ASSERT_TRUE(std::get<0>(BaseLib::getParenthesizedString(
                                    test_string, '{', '}',
                                    test_string.length() - post.length() + 1))
                        .empty());
    }

    // }a{b}"
    for (int i = 0; i < 100; ++i)
    {
        std::string const post = generate_random_alphanumeric_string(20);
        std::string const expected =
            generate_random_alphanumeric_string(20) + "{" + post;
        std::string const test_string = "{" + expected + "}";
        ASSERT_EQ(expected,
                  std::get<0>(BaseLib::getParenthesizedString(test_string, '{',
                                                              '}', 0)));
        ASSERT_TRUE(std::get<0>(BaseLib::getParenthesizedString(
                                    test_string, '{', '}',
                                    test_string.length() - post.length() + 1))
                        .empty());
    }

    // a{b}{c}
    for (int i = 0; i < 100; ++i)
    {
        std::string const pre = generate_random_alphanumeric_string(20);
        std::string const expected_1 = generate_random_alphanumeric_string(20);
        std::string const expected_2 = generate_random_alphanumeric_string(20);
        std::string const test_string =
            pre + "{" + expected_1 + "}{" + expected_2 + "}";
        ASSERT_EQ(expected_1,
                  std::get<0>(BaseLib::getParenthesizedString(test_string, '{',
                                                              '}', 0)));
        ASSERT_EQ(
            expected_2,
            std::get<0>(BaseLib::getParenthesizedString(
                test_string, '{', '}', pre.length() + expected_1.length())));
        ASSERT_TRUE(
            std::get<0>(BaseLib::getParenthesizedString(
                            test_string, '{', '}',
                            test_string.length() - expected_2.length() + 1))
                .empty());
    }
}
