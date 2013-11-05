/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-05-15
 * \brief  Interface tests of global vector classes
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
#include "MathLib/LinAlg/FinalizeVectorAssembly.h"

#ifdef USE_LIS
#include "MathLib/LinAlg/Lis/LisVector.h"
#endif

namespace
{

template <class T_VECTOR>
void checkGlobalVectorInterface()
{
    T_VECTOR x(10);

    ASSERT_EQ(10u, x.size());
    ASSERT_TRUE(x.getRangeBegin()>=0);
    ASSERT_TRUE(x.getRangeEnd()>=0);


    finalizeVectorAssembly(x);

    /*
    ASSERT_EQ(.0, x.get(0));
    x.set(0, 1.0);
    ASSERT_EQ(1.0, x.get(0));
    ASSERT_EQ(0.0, x.get(1));
    x.add(0, 1.0);
    ASSERT_EQ(2.0, x.get(0));

    T_VECTOR y(x);
    ASSERT_EQ(2.0, y.get(0));
    ASSERT_EQ(0.0, y.get(1));
    y += x;
    ASSERT_EQ(4.0, y.get(0));
    y -= x;
    ASSERT_EQ(2.0, y.get(0));
    y = 1.0;
    ASSERT_EQ(1.0, y.get(0));
    y = x;
    ASSERT_EQ(2.0, y.get(0));

    std::vector<double> local_vec(2, 1.0);
    std::vector<std::size_t> vec_pos(2);
    vec_pos[0] = 0;
    vec_pos[1] = 3;
    y.add(vec_pos, local_vec);
    ASSERT_EQ(3.0, y.get(0));
    ASSERT_EQ(0.0, y.get(1));
    ASSERT_EQ(1.0, y.get(3));
    */
}

} // end namespace

TEST(Math, CheckInterface_DenseVector)
{
    checkGlobalVectorInterface<MathLib::DenseVector<double> >();
}

#ifdef USE_LIS
TEST(Math, CheckInterface_LisVector)
{
    checkGlobalVectorInterface<MathLib::LisVector >();
}
#endif

