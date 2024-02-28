/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include "ProcessLib/ConstitutiveRelations/Base.h"

template <class T>
struct ProcessLib_ConstitutiveSettingPrevState_Typed : ::testing::Test
{
};

using ProcessLib_ConstitutiveSettingPrevState_TypedTestCases =
    ::testing::Types<int, std::string>;

TYPED_TEST_SUITE(ProcessLib_ConstitutiveSettingPrevState_Typed,
                 ProcessLib_ConstitutiveSettingPrevState_TypedTestCases);

TYPED_TEST(ProcessLib_ConstitutiveSettingPrevState_Typed, StaticAssertions_Ctor)
{
    using namespace ProcessLib::ConstitutiveRelations;
    using T = TypeParam;

    // copy, move, default
    static_assert(std::is_default_constructible_v<PrevState<T>>);
    static_assert(std::is_copy_constructible_v<PrevState<T>>);
    static_assert(std::is_move_constructible_v<PrevState<T>>);

    // from value type
    static_assert(std::is_constructible_v<PrevState<T>, T>);
    static_assert(std::is_constructible_v<PrevState<T>, T&>);
    static_assert(std::is_constructible_v<PrevState<T>, T const&>);
    static_assert(std::is_constructible_v<PrevState<T>, T&&>);
}

TYPED_TEST(ProcessLib_ConstitutiveSettingPrevState_Typed,
           StaticAssertions_Assignment)
{
    using namespace ProcessLib::ConstitutiveRelations;
    using T = TypeParam;

    // copy, move
    static_assert(std::is_copy_assignable_v<PrevState<T>>);
    static_assert(std::is_move_assignable_v<PrevState<T>>);

    // value type
    static_assert(std::is_assignable_v<PrevState<T>, T>);
    static_assert(std::is_assignable_v<PrevState<T>, T&>);
    static_assert(std::is_assignable_v<PrevState<T>, T const&>);
    static_assert(std::is_assignable_v<PrevState<T>, T&&>);

    // overwriting value type with prev state not possible
    static_assert(!std::is_assignable_v<T, PrevState<T>>);
    static_assert(!std::is_assignable_v<T, PrevState<T>&>);
    static_assert(!std::is_assignable_v<T, PrevState<T> const&>);
    static_assert(!std::is_assignable_v<T, PrevState<T>&&>);
}

TEST(ProcessLib_ConstitutiveSettingPrevState, Tests1)
{
    using P = ProcessLib::ConstitutiveRelations::PrevState<std::string>;

    P p;  // default constructed

    ASSERT_TRUE(p->empty());  // member function access via ->
    ASSERT_EQ("", *p);        // value access via operator*

    *p = "some string";  // assignment to underlying value
    ASSERT_EQ("some string", *p);

    p = "something else";  // direct assignment
    ASSERT_EQ("something else", *p);
}

TEST(ProcessLib_ConstitutiveSettingPrevState, Tests2)
{
    using P = ProcessLib::ConstitutiveRelations::PrevState<std::string>;

    P p{"msg"};

    ASSERT_EQ(3, p->length());
    ASSERT_EQ("msg", *p);

    std::string s = "some string";
    p = s;
    ASSERT_EQ("some string", *p);

    s = "something else";
    p = std::move(s);
    ASSERT_EQ("something else", *p);
}

TEST(ProcessLib_ConstitutiveSettingPrevState, PrevStateOf)
{
    using namespace ProcessLib::ConstitutiveRelations;

    // empty tuple
    {
        using Tuple = std::tuple<>;
        static_assert(std::is_same_v<Tuple, PrevStateOf<Tuple>>);
    }

    // one entry
    {
        using Tuple = std::tuple<int>;
        using TuplePrev = std::tuple<PrevState<int>>;
        static_assert(std::is_same_v<TuplePrev, PrevStateOf<Tuple>>);
    }

    // multiple entries
    {
        using Tuple = std::tuple<int, double, std::string>;
        using TuplePrev = std::
            tuple<PrevState<int>, PrevState<double>, PrevState<std::string>>;
        static_assert(std::is_same_v<TuplePrev, PrevStateOf<Tuple>>);
    }
}

TEST(ProcessLib_ConstitutiveSettingPrevState, PrevStateOfAssign)
{
    using namespace ProcessLib::ConstitutiveRelations;

    // empty tuple
    {
        using Tuple = std::tuple<>;

        Tuple const t;
        PrevStateOf<Tuple> tp;

        assign(tp, t);  // only needs to compile
    }

    // one entry
    {
        using Tuple = std::tuple<int>;
        using TuplePrev = std::tuple<PrevState<int>>;

        Tuple const t{5};
        TuplePrev tp{6};

        assign(tp, t);

        EXPECT_EQ(5, *std::get<0>(tp));
    }

    // multiple entries
    {
        using Tuple = std::tuple<int, double, std::string>;
        using TuplePrev = std::
            tuple<PrevState<int>, PrevState<double>, PrevState<std::string>>;

        Tuple const t{5, 7.5, "hello"};
        TuplePrev tp{6, 8.25, "bye"};

        assign(tp, t);

        EXPECT_EQ(5, *std::get<0>(tp));
        EXPECT_EQ(7.5, *std::get<1>(tp));
        EXPECT_EQ("hello", *std::get<2>(tp));
    }
}
