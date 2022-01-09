/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include <tuple>

#include "BaseLib/TMP.h"

TEST(BaseLibTMP, Concat)
{
    namespace T = BaseLib::TMP;

    // empty lists
    {
        using Expected = std::tuple<>;
        using First = std::tuple<>;
        using Second = std::tuple<>;
        using Actual = T::Concat_t<First, Second>;
        static_assert(std::is_same_v<Expected, Actual>);
    }

    // single-element and empty
    {
        using Expected = std::tuple<double>;
        using First = std::tuple<double>;
        using Second = std::tuple<>;
        using Actual = T::Concat_t<First, Second>;
        static_assert(std::is_same_v<Expected, Actual>);
    }

    // empty and single-element
    {
        using Expected = std::tuple<int>;
        using First = std::tuple<>;
        using Second = std::tuple<int>;
        using Actual = T::Concat_t<First, Second>;
        static_assert(std::is_same_v<Expected, Actual>);
    }

    // twice the same
    {
        using Expected = std::tuple<unsigned, unsigned>;
        using First = std::tuple<unsigned>;
        using Second = std::tuple<unsigned>;
        using Actual = T::Concat_t<First, Second>;
        static_assert(std::is_same_v<Expected, Actual>);
    }

    // two lists
    {
        using Expected = std::tuple<double, std::string, int, char>;
        using First = std::tuple<double, std::string>;
        using Second = std::tuple<int, char>;
        using Actual = T::Concat_t<First, Second>;
        static_assert(std::is_same_v<Expected, Actual>);
    }
}

TEST(BaseLibTMP, Filter)
{
    namespace T = BaseLib::TMP;

    // empty list
    {
        using Expected = T::List<>;
        using List = std::tuple<>;
        auto pred = [](auto*) { return true; };
        using Actual = decltype(T::filter<List>(pred));
        static_assert(std::is_same_v<Expected, Actual>);
    }

    // keep all
    {
        using Expected = T::List<int, double, std::string>;
        using List = std::tuple<int, double, std::string>;
        auto pred = [](auto*) { return true; };
        using Actual = decltype(T::filter<List>(pred));
        static_assert(std::is_same_v<Expected, Actual>);
    }

    // keep none
    {
        using Expected = T::List<>;
        using List = std::tuple<int, double, std::string>;
        auto pred = [](auto*) { return false; };
        using Actual = decltype(T::filter<List>(pred));
        static_assert(std::is_same_v<Expected, Actual>);
    }

    // keep some
    {
        using Expected = T::List<int, unsigned, char>;
        using List = std::tuple<int, double, unsigned, std::string, char>;
        auto pred = []<typename T>(T*) { return std::is_integral_v<T>; };
        using Actual = decltype(T::filter<List>(pred));
        static_assert(std::is_same_v<Expected, Actual>);
    }
}

TEST(BaseLibTMP, Map)
{
    namespace T = BaseLib::TMP;

    // empty list
    {
        using Expected = std::tuple<>;
        using List = std::tuple<>;
        using Actual = T::Map_t<std::add_pointer_t, List>;
        static_assert(std::is_same_v<Expected, Actual>);
    }

    // non-empty list
    {
        using Expected = std::tuple<int*, double*, std::string*>;
        using List = std::tuple<int, double, std::string>;
        using Actual = T::Map_t<std::add_pointer_t, List>;
        static_assert(std::is_same_v<Expected, Actual>);
    }
}

TEST(BaseLibTMP, MapToArray)
{
    namespace T = BaseLib::TMP;

    // empty list
    {
        using Expected = std::array<std::nullptr_t, 0>;
        using List = std::tuple<>;
        auto mapping = [](auto*) { return 0; };
        auto actual = T::map_to_array<List>(mapping);
        static_assert(std::is_same_v<Expected, decltype(actual)>);
    }

    // non-empty list
    {
        std::array<bool, 5> expected{true, false, true, false, true};
        using List = std::tuple<int, double, unsigned, std::string, char>;
        auto mapping = []<typename T>(T*) { return std::is_integral_v<T>; };
        auto actual = T::map_to_array<List>(mapping);
        EXPECT_THAT(actual, testing::ContainerEq(expected));
    }
}

TEST(BaseLibTMP, ForEach)
{
    namespace T = BaseLib::TMP;

    // empty list
    {
        std::vector<bool> expected{};

        using List = std::tuple<>;
        std::vector<bool> actual{};
        auto callback = [&actual](auto*) { actual.push_back(true); };

        T::foreach<List>(callback);

        EXPECT_TRUE(actual.empty());
    }

    // non-empty list
    {
        std::vector<bool> expected{true, false, true, false, true};

        using List = std::tuple<int, double, unsigned, std::string, char>;
        std::vector<bool> actual{};
        auto callback = [&actual]<typename T>(T*)
        { actual.push_back(std::is_integral_v<T>); };

        T::foreach<List>(callback);

        EXPECT_THAT(actual, testing::ContainerEq(expected));
    }
}

TEST(BaseLibTMP, Contains)
{
    namespace T = BaseLib::TMP;

    // empty list
    {
        using List = std::tuple<>;
        using Type = double;
        auto constexpr actual = T::contains<List, Type>();
        static_assert(!actual);
    }

    // found
    {
        using List = std::tuple<int, double, unsigned, std::string, char>;
        static_assert(T::contains<List, int>());
        static_assert(T::contains<List, double>());
        static_assert(T::contains<List, char>());
    }

    // not found
    {
        using List = std::tuple<int, double, unsigned, std::string, char>;
        static_assert(!T::contains<List, unsigned char>());
    }
}
