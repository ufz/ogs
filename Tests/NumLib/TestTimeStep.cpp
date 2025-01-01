/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include "NumLib/TimeStepping/TimeStep.h"

TEST(NumLib, TimeStep)
{
    // initial
    NumLib::TimeStep t1{NumLib::Time{1.}};
    ASSERT_EQ(1., t1.current()());
    ASSERT_EQ(1., t1.previous()());
    ASSERT_EQ(.0, t1.dt());
    ASSERT_EQ(0u, t1.timeStepNumber());

    NumLib::TimeStep t2{NumLib::Time{0.}, NumLib::Time{1.}, 1};
    ASSERT_EQ(1., t2.current()());
    ASSERT_EQ(0., t2.previous()());
    ASSERT_EQ(1., t2.dt());
    ASSERT_EQ(1u, t2.timeStepNumber());

    // copy
    const NumLib::TimeStep& t3(t2);
    ASSERT_EQ(1., t3.current()());
    ASSERT_EQ(0., t3.previous()());
    ASSERT_EQ(1., t3.dt());
    ASSERT_EQ(1u, t3.timeStepNumber());

    // comparison
    ASSERT_EQ(t2.current(), t3.current());
}
