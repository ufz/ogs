/**
 * \file
 * \author Norihiro Watanabe
 * \author Wenqing Wang
 * \date   2013-05-15 -- 2014-01
 * \brief  Interface tests of global vector classes
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
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

#ifdef USE_PETSC
#include "MathLib/LinAlg/PETSc/PETScVector.h"
#endif

namespace
{

template <class T_VECTOR>
void checkGlobalVectorInterface()
{
    T_VECTOR x(10);

    ASSERT_EQ(10u, x.size());
    ASSERT_EQ(0u, x.getRangeBegin());
    ASSERT_EQ(10u, x.getRangeEnd());

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

#ifdef USE_PETSC // or MPI
template <class T_VECTOR>
void checkGlobalVectorInterfaceMPI()
{
    int msize;
    MPI_Comm_size(PETSC_COMM_WORLD, &msize);

    ASSERT_EQ(3u, msize);

    // -------------------------------------------------------------------
    // PETSc determined partitioning
    T_VECTOR x(16);

    ASSERT_EQ(16u, x.size());
    ASSERT_EQ(x.getRangeEnd()-x.getRangeBegin(), x.getLocalSize());

    const int r0 = x.getRangeBegin();
    //x.get(0) is expensive, only get local value. Use it for test purpose
    ASSERT_EQ(.0, x.get(r0));

    x = 10.;

    // Value of x is not copied to y
    const bool deep_copy = false;
    T_VECTOR y(x, deep_copy);
    ASSERT_EQ(0, y.get(r0));

    y = 10.0;
    ASSERT_EQ(10, y.get(r0));

    y += x;
    ASSERT_EQ(20, y.get(r0));
    ASSERT_EQ(80., y.getNorm());

    y -= x;
    ASSERT_EQ(10, y.get(r0));
    ASSERT_EQ(40., y.getNorm());

    std::vector<double> local_vec(2, 10.0);
    std::vector<int> vec_pos(2);

    vec_pos[0] = r0;   // any index in [0,15]
    vec_pos[1] = r0+1; // any index in [0,15]

    y.add(vec_pos, local_vec);

    double normy = std::sqrt(6.0*400+10.0*100);

    ASSERT_NEAR(0.0, normy-y.getNorm(), 1.e-10);

    double x0[16];
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

    y.getGlobalVector(x0);

    ASSERT_ARRAY_NEAR(z, x0, 16, 1e-10);

    // -------------------------------------------------------------------
    // User determined partitioning
    const bool is_gloabal_size = false;
    T_VECTOR x_fixed_p(2, is_gloabal_size);

    ASSERT_EQ(6u, x_fixed_p.size());

    int mrank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &mrank);

    ASSERT_EQ(2*mrank, x_fixed_p.getRangeBegin());
    ASSERT_EQ(2*mrank+2, x_fixed_p.getRangeEnd());

    vec_pos[0] = 2 * mrank;
    vec_pos[1] = vec_pos[0] + 1;
    local_vec[0] = 1.;
    local_vec[1] = 2.;
    for(unsigned i=0; i<3; i++)
    {
        const unsigned j = 2 * i;
        z[j] = 1.0;
        z[j+1] = 2.0;
    }
    x_fixed_p.set(vec_pos, local_vec);
    x_fixed_p.getGlobalVector(x0);

    ASSERT_ARRAY_NEAR(z, x0, 6, 1e-10);

    // check local array
    double *loc_v = x_fixed_p.getLocalVector();
    z[0] = 1.0;
    z[1] = 2.0;

    ASSERT_ARRAY_NEAR(z, loc_v, 2, 1e-10);

    // Deep copy
    MathLib::finalizeVectorAssembly(x_fixed_p);
    T_VECTOR x_deep_copied(x_fixed_p);
    ASSERT_NEAR(sqrt(3.0*5), x_deep_copied.getNorm(), 1.e-10);
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
#ifdef USE_PETSC
TEST(MPITest_Math, CheckInterface_PETScVector)
{
    checkGlobalVectorInterfaceMPI<MathLib::PETScVector >();
}
#endif
