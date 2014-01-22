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
#include "../TestTools.h"

#include "MathLib/LinAlg/Dense/DenseVector.h"
#include "MathLib/LinAlg/FinalizeVectorAssembly.h"

#ifdef USE_LIS
#include "MathLib/LinAlg/Lis/LisVector.h"
#endif

#ifdef OGS_USE_PETSC
#include "MathLib/LinAlg/PETSc/PETScVector.h"
#include "MathLib/LinAlg/PETSc/InfoMPI.h"
#endif

namespace
{

template <class T_VECTOR>
void checkGlobalVectorInterface()
{
    T_VECTOR x(10);

    ASSERT_EQ(10u, x.size());

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
}

#ifdef OGS_USE_PETSC // or MPI
template <class T_VECTOR>
void checkGlobalVectorInterfaceMPI()
{
    ASSERT_EQ(3u, BaseLib:: InfoMPI::getSize());

    T_VECTOR x(16);

    const int r0 = x.getRangeBegin();

    ASSERT_EQ(16u, x.size());
    ASSERT_TRUE(r0 >= 0);
    ASSERT_TRUE(x.getRangeEnd() >= 0);

    //x.get(0) is expensive, only get local value. Use it for test purpose
    ASSERT_EQ(.0, x.get(r0));

    x = 10.;

    // Value of x is not copied to y
    T_VECTOR y(x);

    y = 10.0;
    y += x;
    //y.get(0) is expensive
    ASSERT_EQ(20, y.get(r0));
    ASSERT_EQ(80., y.getNorm());

    y -= x;
    //y.get(0) is expensive
    ASSERT_EQ(10, y.get(r0));
    ASSERT_EQ(40., y.getNorm());

    std::vector<double> local_vec(2, 10.0);
    std::vector<int> vec_pos(2);

    vec_pos[0] = r0;   // any index in [0,15]
    vec_pos[1] = r0+1; // any index in [0,15]

    y.add(vec_pos, local_vec);

    double normy = std::sqrt(6.0*400+10.0*100);

    ASSERT_NEAR(normy-y.getNorm(), 0.0, 1.e-10);

    double x0[16];
    double x1[16];
    double z[] =
    {
        2.0000000000000000e+01,
        2.0000000000000000e+01,
        1.0000000000000000e+01,
        1.0000000000000000e+01,
        1.0000000000000000e+01,
        1.0000000000000000e+01,
        2.0000000000000000e+01,
        2.0000000000000000e+01,
        1.0000000000000000e+01,
        1.0000000000000000e+01,
        1.0000000000000000e+01,
        2.0000000000000000e+01,
        2.0000000000000000e+01,
        1.0000000000000000e+01,
        1.0000000000000000e+01,
        1.0000000000000000e+01
    };

    y.getGlobalEntries(x0, x1);

    ASSERT_ARRAY_NEAR(x0, z, 16, 1e-10);
    ASSERT_ARRAY_NEAR(x1, z, 16, 1e-10);
}
#endif

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

//--------------------------------------------
#ifdef OGS_USE_PETSC
TEST(Math, CheckInterface_PETScVector)
{
    checkGlobalVectorInterfaceMPI<MathLib::PETScVector >();
}
#endif


