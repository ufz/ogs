/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-05-15
 * \brief  Implementation tests of Vector classes.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include "MathLib/LinAlg/Dense/DenseVector.h"

TEST(Math, CheckInterface_DenseVector)
{
    MathLib::DenseVector<double> x(10);

    ASSERT_EQ(10u, x.size());
    ASSERT_EQ(0u, x.getRangeBegin());
    ASSERT_EQ(10u, x.getRangeEnd());

    ASSERT_EQ(.0, x.get(0));
    x.set(0, 1.0);
    ASSERT_EQ(1.0, x.get(0));
    x.add(0, 1.0);
    ASSERT_EQ(2.0, x.get(0));

    MathLib::DenseVector<double> y(x);
    ASSERT_EQ(2.0, y.get(0));
    y += x;
    ASSERT_EQ(4.0, y.get(0));
    y -= x;
    ASSERT_EQ(2.0, y.get(0));
    y = 1.0;
    ASSERT_EQ(1.0, y.get(0));
    y = x;
    ASSERT_EQ(2.0, y.get(0));

}


