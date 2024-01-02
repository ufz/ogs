/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ProcessLib/Utils/TransposeInPlace.h"

TEST(ProcessLib, TransposeInplaceFixed)
{
    std::vector<double> const values_in{
        1, 2,  3,  4,  //
        5, 6,  7,  8,  //
        9, 10, 11, 12  //
    };
    constexpr unsigned num_rows = 3;

    std::vector<double> const values_out_expected{
        1, 5, 9,   //
        2, 6, 10,  //
        3, 7, 11,  //
        4, 8, 12   //
    };

    std::vector<double> values_out_actual(values_in);

    ProcessLib::transposeInPlace<num_rows>(values_out_actual);

    ASSERT_THAT(values_out_actual, testing::ContainerEq(values_out_expected));
}

TEST(ProcessLib, TransposeInplaceDynamic)
{
    std::vector<double> const values_in{
        1, 2,  3,  4,  //
        5, 6,  7,  8,  //
        9, 10, 11, 12  //
    };
    unsigned const num_rows = 3;

    std::vector<double> const values_out_expected{
        1, 5, 9,   //
        2, 6, 10,  //
        3, 7, 11,  //
        4, 8, 12   //
    };

    std::vector<double> values_out_actual(values_in);

    ProcessLib::transposeInPlace(values_out_actual, num_rows);

    ASSERT_THAT(values_out_actual, testing::ContainerEq(values_out_expected));
}
