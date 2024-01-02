/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include <array>

#include "BaseLib/StrongType.h"

TEST(BaseLibStrongType, CopyOrMoveValueIn_Class)
{
    using StrongVector =
        BaseLib::StrongType<std::vector<int>, struct StrongVectorTag>;

    std::vector<int> const v{1, 2};

    // copy value in
    {
        StrongVector const sv{v};

        EXPECT_THAT(sv(), ::testing::ContainerEq(v));
    }

    // copy value in, mutable vector
    {
        std::vector<int> v_copy = v;

        StrongVector sv{v_copy};

        EXPECT_THAT(sv(), ::testing::ContainerEq(v));

        sv()[0] = 8;
        sv().back() = 6;

        EXPECT_THAT(v_copy, ::testing::ContainerEq(v))
            << "the input vector was modified";
        EXPECT_EQ(8, sv().front());
        EXPECT_EQ(6, sv()[1]);
    }

    // move value in
    {
        using StrongIntPtr =
            BaseLib::StrongType<std::shared_ptr<int>, struct StrongIntPtrTag>;

        auto ip = std::make_shared<int>(23);
        StrongIntPtr const sip{std::move(ip)};

        EXPECT_FALSE(ip) << "the input pointer was not moved from";
        EXPECT_EQ(23, **sip);
    }
}

TEST(BaseLibStrongType, CopyValueIn_PrimitiveType)
{
    using StrongDouble = BaseLib::StrongType<double, struct StrongDoubleTag>;

    double const d = 3;

    // copy value in
    {
        StrongDouble sd{d};

        EXPECT_EQ(3, sd());
    }

    // copy value in, mutable value
    {
        double d2 = d;

        StrongDouble sd{d2};

        EXPECT_EQ(3, sd());

        sd() = 5;

        EXPECT_EQ(3, d2) << "the original data was modified";
        EXPECT_EQ(5, sd());
    }
}

TEST(BaseLibStrongType, Construction)
{
    using Vector = std::vector<int>;
    using StrongVector = BaseLib::StrongType<Vector, struct StrongVectorTag>;

    // default construct
    {
        StrongVector const sv;

        EXPECT_TRUE(sv().empty());
    }

    StrongVector const sv{{1, 2}};

    // copy construct
    {
        StrongVector sv2{sv};
        StrongVector sv3{sv2};

        EXPECT_THAT(sv2(), ::testing::ContainerEq(sv()));
        EXPECT_THAT(sv3(), ::testing::ContainerEq(sv()));

        sv3()[0] = 8;
        sv3().back() = 6;

        EXPECT_THAT(sv2(), ::testing::ContainerEq(sv()))
            << "only sv3 should have changed, but sv2 has, too";
    }

    // move construct
    {
        using StrongIntPtr =
            BaseLib::StrongType<std::shared_ptr<int>, struct StrongIntPtrTag>;

        StrongIntPtr sip{std::make_shared<int>(42)};
        StrongIntPtr const sip2{std::move(sip)};

        EXPECT_FALSE(*sip);
        EXPECT_EQ(42, **sip2);
    }
}

TEST(BaseLibStrongType, Assignment)
{
    using StrongVector =
        BaseLib::StrongType<std::vector<int>, struct StrongVectorTag>;

    StrongVector const sv{{7, 9}};

    // copy assign
    {
        StrongVector sv2{sv};
        StrongVector sv3;
        sv3 = sv2;

        EXPECT_THAT(sv3(), ::testing::ContainerEq(sv()));

        sv3()[0] = 8;
        sv3().back() = 6;

        EXPECT_THAT(sv2(), ::testing::ContainerEq(sv()))
            << "only sv3 should have changed, but sv2 has, too";
    }

    // move assign
    {
        using StrongIntPtr =
            BaseLib::StrongType<std::shared_ptr<int>, struct StrongIntPtrTag>;

        StrongIntPtr sip{std::make_shared<int>(42)};
        StrongIntPtr sip2;
        sip2 = std::move(sip);

        EXPECT_FALSE(*sip);
        EXPECT_EQ(42, **sip2);
    }
}

TEST(BaseLibStrongType, AccessOperators)
{
    using T = std::vector<int>;
    using ST = BaseLib::StrongType<T, struct Tag>;

    // const access
    {
        ST const st{{1, 2, 3, 4}};

        EXPECT_EQ(4, st().size());
        EXPECT_EQ(4, (*st).size());
        EXPECT_EQ(4, st->size());

        EXPECT_THAT(st(), ::testing::ContainerEq(*st));
        EXPECT_EQ(st(), *st);  // uses the std::vector comparison operator
    }

    // non-const access
    {
        T const t{1, 2, 3, 4};
        ST st{t};

        EXPECT_EQ(4, st->size());

        st().push_back(7);
        st().push_back(0);

        EXPECT_EQ(6, (*st).size());

        (*st).pop_back();

        EXPECT_EQ(5, st->size());

        st->pop_back();

        EXPECT_EQ(4, st().size());
        EXPECT_THAT(*st, ::testing::ContainerEq(t));
    }
}

TEST(BaseLibStrongType, Constexpr)
{
    using ST = BaseLib::StrongType<int, struct Tag>;
    constexpr ST s{5};

    static_assert(s() == 5);

    std::array<double, s()> a = {};
    static_assert(a.size() == 5);
}

TEST(BaseLibStrongType, CompileTimePure)
{
    auto check = []<typename T>(std::type_identity<T>)
    {
        using ST = BaseLib::StrongType<T, struct Tag1>;

        // construction and assignment
        {
            static_assert(std::is_default_constructible_v<ST> ==
                          std::is_default_constructible_v<T>);
            static_assert(std::is_copy_constructible_v<ST> ==
                          std::is_copy_constructible_v<T>);
            static_assert(std::is_move_constructible_v<ST> ==
                          std::is_move_constructible_v<T>);

            static_assert(std::is_copy_assignable_v<ST> ==
                          std::is_copy_assignable_v<T>);
            static_assert(std::is_move_assignable_v<ST> ==
                          std::is_move_assignable_v<T>);
        }

        // construction from value type
        {
            static_assert(std::is_constructible_v<ST, T>);
            static_assert(std::is_constructible_v<ST, T&>);
            static_assert(std::is_constructible_v<ST, T const&>);
            static_assert(std::is_constructible_v<ST, T&&>);
        }

        // no implicit conversions to value type
        {
            static_assert(!std::is_convertible_v<ST, T>);
            static_assert(!std::is_convertible_v<ST, T&>);
            static_assert(!std::is_convertible_v<ST, T const&>);
            static_assert(!std::is_convertible_v<ST, T&&>);
        }

        // self test of the unit test: normal function calls work
        {
            static_assert(std::is_invocable_v<void (*)(ST), ST> ==
                          std::is_move_constructible_v<ST>);
            // or, equivalently (passing r-value ref instead of value):
            static_assert(std::is_invocable_v<void (*)(ST), ST&&> ==
                          std::is_move_constructible_v<ST>);

            static_assert(std::is_invocable_v<void (*)(ST), ST const&> ==
                          std::is_copy_constructible_v<ST>);

            static_assert(std::is_invocable_v<void (*)(ST&), ST&>);
            static_assert(std::is_invocable_v<void (*)(ST const&), ST>);
            static_assert(std::is_invocable_v<void (*)(ST &&), ST>);
        }
    };

    check(std::type_identity<std::vector<int>>{});
    check(std::type_identity<std::unique_ptr<char>>{});

    struct NoMoveCtor
    {
        NoMoveCtor() = default;
        NoMoveCtor(NoMoveCtor&&) = delete;
    };

    check(std::type_identity<NoMoveCtor>{});
}

TEST(BaseLibStrongType, CompileTimeMixed)
{
    using V1 = std::vector<int>;
    using S1 = BaseLib::StrongType<V1, struct Tag1>;
    using S2 = BaseLib::StrongType<V1, struct Tag2>;

    // no mixed construction possible
    {
        static_assert(!std::is_constructible_v<S1, S2>);
        static_assert(!std::is_constructible_v<S1, S2&>);
        static_assert(!std::is_constructible_v<S1, S2 const&>);
        static_assert(!std::is_constructible_v<S1, S2&&>);
    }

    // no mixed implicit conversions
    {
        static_assert(!std::is_convertible_v<S1, S2>);
        static_assert(!std::is_convertible_v<S1, S2&>);
        static_assert(!std::is_convertible_v<S1, S2 const&>);
        static_assert(!std::is_convertible_v<S1, S2&&>);
    }

    // no mixed assignment
    {
        static_assert(!std::is_assignable_v<S1&, S2>);
        static_assert(!std::is_assignable_v<S1&, S2&>);
        static_assert(!std::is_assignable_v<S1&, S2 const&>);
        static_assert(!std::is_assignable_v<S1&, S2&&>);
    }

    // no mixed function calls (this is a central property of strong types)
    {
        static_assert(!std::is_invocable_v<void (*)(S2), S1>);
        static_assert(!std::is_invocable_v<void (*)(S2&), S1>);
        static_assert(!std::is_invocable_v<void (*)(S2 const&), S1>);
        static_assert(!std::is_invocable_v<void (*)(S2 &&), S1>);

        static_assert(!std::is_invocable_v<void (*)(S2), S1&>);
        static_assert(!std::is_invocable_v<void (*)(S2&), S1&>);
        static_assert(!std::is_invocable_v<void (*)(S2 const&), S1&>);
        static_assert(!std::is_invocable_v<void (*)(S2 &&), S1&>);

        static_assert(!std::is_invocable_v<void (*)(S2), S1 const&>);
        static_assert(!std::is_invocable_v<void (*)(S2&), S1 const&>);
        static_assert(!std::is_invocable_v<void (*)(S2 const&), S1 const&>);
        static_assert(!std::is_invocable_v<void (*)(S2 &&), S1 const&>);

        static_assert(!std::is_invocable_v<void (*)(S2), S1&&>);
        static_assert(!std::is_invocable_v<void (*)(S2&), S1&&>);
        static_assert(!std::is_invocable_v<void (*)(S2 const&), S1&&>);
        static_assert(!std::is_invocable_v<void (*)(S2 &&), S1&&>);
    }

    // no function calls with the base type (this is a central property of
    // strong types)
    {
        static_assert(!std::is_invocable_v<void (*)(S2), V1>);
        static_assert(!std::is_invocable_v<void (*)(S2&), V1>);
        static_assert(!std::is_invocable_v<void (*)(S2 const&), V1>);
        static_assert(!std::is_invocable_v<void (*)(S2 &&), V1>);

        static_assert(!std::is_invocable_v<void (*)(S2), V1&>);
        static_assert(!std::is_invocable_v<void (*)(S2&), V1&>);
        static_assert(!std::is_invocable_v<void (*)(S2 const&), V1&>);
        static_assert(!std::is_invocable_v<void (*)(S2 &&), V1&>);

        static_assert(!std::is_invocable_v<void (*)(S2), V1 const&>);
        static_assert(!std::is_invocable_v<void (*)(S2&), V1 const&>);
        static_assert(!std::is_invocable_v<void (*)(S2 const&), V1 const&>);
        static_assert(!std::is_invocable_v<void (*)(S2 &&), V1 const&>);

        static_assert(!std::is_invocable_v<void (*)(S2), V1&&>);
        static_assert(!std::is_invocable_v<void (*)(S2&), V1&&>);
        static_assert(!std::is_invocable_v<void (*)(S2 const&), V1&&>);
        static_assert(!std::is_invocable_v<void (*)(S2 &&), V1&&>);
    }
}

TEST(BaseLibStrongType, DefaultInitialize)
{
    // int is a primitive type, but it should be initialized properly, anyway
    BaseLib::StrongType<int, struct Tag> si;

    ASSERT_EQ(0, si());
}
