/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <vector>
#include <functional>
#include <numeric>
#include "NumLib/Assembler/SerialExecutor.h"

template <typename ContainerElement_>
class NumLibSerialExecutor : public ::testing::Test
{
public:
    using ContainerElement = ContainerElement_;
    using Container = std::vector<ContainerElement>;
    using PtrContainer = std::vector<ContainerElement*>;

    template<typename Callback>
    void test(Callback const& cb)
    {
        Container reference(size);
        std::iota(reference.begin(), reference.end(), 0);

        Container container_back(reference);

        PtrContainer container;
        container.reserve(size);
        for (auto& el : container_back) container.push_back(&el);

        cb(container, reference);

        ASSERT_TRUE(referenceIsZero(reference));
    }

    static void subtractFromReferenceStatic(
            ContainerElement const value, std::size_t const index,
            Container& reference)
    {
        reference[index] -= value;
    }

    void subtractFromReference(
            ContainerElement const value, std::size_t const index,
            Container& reference) const
    {
        reference[index] -= value;
    }

    static bool
    referenceIsZero(Container const& reference)
    {
        return std::all_of(reference.begin(), reference.end(),
            [](ContainerElement const reference_value)
            {
                return reference_value == 0;
            });
    }

    static std::size_t const size = 100;
};

typedef ::testing::Types<int> TestCases;

TYPED_TEST_CASE(NumLibSerialExecutor, TestCases);

TYPED_TEST(NumLibSerialExecutor, ContainerArgument)
{
    using Elem         = typename TestFixture::ContainerElement;
    using Container    = typename TestFixture::Container;
    using PtrContainer = typename TestFixture::PtrContainer;

    TestFixture::test(
        [](PtrContainer const& ctnr, Container& ref) {
            auto cb_static =
                [](Elem const value, std::size_t const index, Container& ref_inner) {
                    TestFixture::subtractFromReferenceStatic(value, index, ref_inner);
                };

            NumLib::SerialExecutor::executeDereferenced(
                cb_static, ctnr, ref);
        }
    );

    TestFixture::test(
        [this](PtrContainer const& ctnr, Container& ref) {
            NumLib::SerialExecutor::executeMemberDereferenced(
                *static_cast<TestFixture*>(this), &TestFixture::subtractFromReference,
                ctnr, ref
            );
        }
    );
}
