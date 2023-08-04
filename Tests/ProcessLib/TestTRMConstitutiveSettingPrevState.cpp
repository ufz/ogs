/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Base.h"

template <class T>
struct ProcessLib_TRMConstitutiveSettingPrevState_Typed : ::testing::Test
{
};

using ProcessLib_TRMConstitutiveSettingPrevState_TypedTestCases =
    ::testing::Types<int, std::string>;

TYPED_TEST_SUITE(ProcessLib_TRMConstitutiveSettingPrevState_Typed,
                 ProcessLib_TRMConstitutiveSettingPrevState_TypedTestCases);

TYPED_TEST(ProcessLib_TRMConstitutiveSettingPrevState_Typed,
           StaticAssertions_Ctor)
{
    using namespace ProcessLib::ThermoRichardsMechanics;
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

TYPED_TEST(ProcessLib_TRMConstitutiveSettingPrevState_Typed,
           StaticAssertions_Assignment)
{
    using namespace ProcessLib::ThermoRichardsMechanics;
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

TEST(ProcessLib_TRMConstitutiveSettingPrevState, Tests1)
{
    using P = ProcessLib::ThermoRichardsMechanics::PrevState<std::string>;

    P p;  // default constructed

    ASSERT_TRUE(p->empty());  // member function access via ->
    ASSERT_EQ("", *p);        // value access via operator*

    *p = "some string";  // assignment to underlying value
    ASSERT_EQ("some string", *p);

    p = "something else";  // direct assignment
    ASSERT_EQ("something else", *p);
}

TEST(ProcessLib_TRMConstitutiveSettingPrevState, Tests2)
{
    using P = ProcessLib::ThermoRichardsMechanics::PrevState<std::string>;

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
