/**
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <vector>
#include <functional>
#include "AssemblerLib/SerialExecutor.h"

template <typename Container>
class AssemblerLibSerialExecutor : public ::testing::Test
{
    public:
    AssemblerLibSerialExecutor()
    {
        reference.resize(size);
        std::iota(reference.begin(), reference.end(), 0);

        std::copy(reference.begin(), reference.end(),
            std::back_inserter(container));
    }

    void
    subtractFromReference(int const value, std::size_t const index)
    {
        reference[index] -= value;
    }

    bool
    referenceIsZero() const
    {
        return std::all_of(reference.begin(), reference.end(),
            [](int const reference_value)
            {
                return reference_value == 0;
            });
    }

    static std::size_t const size = 100;
    std::vector<int> reference;
    Container container;
};

TYPED_TEST_CASE_P(AssemblerLibSerialExecutor);

using namespace std::placeholders;

TYPED_TEST_P(AssemblerLibSerialExecutor, ContainerArgument)
{
    AssemblerLib::SerialExecutor::execute(
        std::bind(&TestFixture::subtractFromReference, this, _1, _2),
        this->container);
    ASSERT_TRUE(this->referenceIsZero());
}

REGISTER_TYPED_TEST_CASE_P(AssemblerLibSerialExecutor,
    ContainerArgument);


typedef ::testing::Types <
      std::vector<int>
    > TestTypes;

INSTANTIATE_TYPED_TEST_CASE_P(templated, AssemblerLibSerialExecutor,
    TestTypes);
