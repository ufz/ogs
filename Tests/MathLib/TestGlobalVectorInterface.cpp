/**
 * \file
 * \author Norihiro Watanabe
 * \author Wenqing Wang
 * \date   2013-05-15 -- 2014-01
 * \brief  Interface tests of global vector classes
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>
#include "Tests/TestTools.h"

#include "MathLib/LinAlg/LinAlg.h"

#if defined(USE_PETSC)
#include "MathLib/LinAlg/PETSc/PETScVector.h"
#elif defined(OGS_USE_EIGEN)
#include "MathLib/LinAlg/Eigen/EigenVector.h"
#endif

#include "MathLib/LinAlg/FinalizeVectorAssembly.h"
#include "NumLib/NumericsConfig.h"

using namespace MathLib::LinAlg;

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
    // y += x;
    axpy(y, 1., x);

    ASSERT_EQ(4.0, y.get(0));
    //y -= x;
    axpy(y, -1., x);
    ASSERT_EQ(2.0, y.get(0));
    //y = 1.0;
    set(y, 1.0);
    ASSERT_EQ(1.0, y.get(0));
    // y = x;
    copy(x, y);
    ASSERT_EQ(2.0, y.get(0));

    std::vector<double> local_vec(2, 1.0);
    std::vector<GlobalIndexType> vec_pos(2);
    vec_pos[0] = 0;
    vec_pos[1] = 3;
    y.add(vec_pos, local_vec);
    ASSERT_EQ(3.0, y.get(0));
    ASSERT_EQ(0.0, y.get(1));
    ASSERT_EQ(1.0, y.get(3));
}

#ifdef USE_PETSC
template <class T_VECTOR>
void checkGlobalVectorInterfacePETSc()
{
    int msize;
    MPI_Comm_size(PETSC_COMM_WORLD, &msize);

    ASSERT_EQ(3u, msize);

    // -----------------------------------------------------------------
    // PETSc determined partitioning
    T_VECTOR x(16);

    ASSERT_EQ(16u, x.size());
    ASSERT_EQ(x.getRangeEnd()-x.getRangeBegin(), x.getLocalSize());

    const int r0 = x.getRangeBegin();
    //x.get(0) is expensive, only get local value. Use it for test purpose
    ASSERT_EQ(.0, x.get(r0));

    set(x, 10.);


    // Value of x is not copied to y
    const bool deep_copy = false;
    T_VECTOR y(x, deep_copy);
    ASSERT_EQ(0, y.get(r0));

    set(y, 10.0);
    ASSERT_EQ(10, y.get(r0));

    // y += x
    axpy(y, 1., x);
    ASSERT_EQ(20, y.get(r0));
    ASSERT_EQ(80., norm2(y));

    // y -= x
    axpy(y, -1., x);
    ASSERT_EQ(10, y.get(r0));
    ASSERT_EQ(40., norm2(y));

    std::vector<double> local_vec(2, 10.0);
    std::vector<GlobalIndexType> vec_pos(2);

    vec_pos[0] = r0;   // any index in [0,15]
    vec_pos[1] = r0+1; // any index in [0,15]

    y.add(vec_pos, local_vec);

    double normy = std::sqrt(6.0*400+10.0*100);

    ASSERT_NEAR(0.0, normy-norm2(y), 1.e-10);

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

    // -----------------------------------------------------------------
    // User determined partitioning
    const bool is_global_size = false;
    T_VECTOR x_fixed_p(2, is_global_size);

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
    std::vector<double> loc_v(  x_fixed_p.getLocalSize()
                              + x_fixed_p.getGhostSize() );
    x_fixed_p.copyValues(loc_v);
    z[0] = 1.0;
    z[1] = 2.0;

    ASSERT_ARRAY_NEAR(z, loc_v, 2, 1e-10);

    // Deep copy
    MathLib::finalizeVectorAssembly(x_fixed_p);
    T_VECTOR x_deep_copied(x_fixed_p);
    ASSERT_NEAR(sqrt(3.0*5), norm2(x_deep_copied), 1.e-10);

    // -----------------------------------------------------------------
    // Vector with ghost entries
    /*
         Assume there is a vector distributed over three processes as
          -- rank0 --    --- rank1 ---   -- rank2 --
           0  1  2  3    4  5  6  7  8   9   10   11 
         where the numbers are the global entry indices.
         In each trunk of entries of a rank, there are ghost entries in
         other ranks attached and their global entry indices are:
         rank0: 6 8 10
         rank1: 0 9
         rank2: 3 5

         Assuming the values of the entries are just their global indices,
         we have local arrays as:
         rank0: 0 1 2 3     6 8 10
         rank1: 4 5 6 7 8   0 9
         rank2: 9 10 11     3 5

         The above ghost entry embedded vector is realized by the following
         test.
    */    
    std::size_t local_vec_size = 4;
    if (mrank == 1)
        local_vec_size = 5;
    else if (mrank == 2)
        local_vec_size = 3;
    std::vector<GlobalIndexType> non_ghost_ids(local_vec_size);
    std::vector<double> non_ghost_vals(local_vec_size);
    std::size_t nghosts = 3;
    if (mrank)
        nghosts = 2;    
    std::vector<GlobalIndexType> ghost_ids(nghosts);
    std::vector<double> expected;
    switch (mrank)
    {
        case 0:
            non_ghost_ids  = {0, 1, 2, 3};
            non_ghost_vals = {0., 1., 2., 3.};
            ghost_ids      = {6, 8, 10};
            expected       = {0., 1., 2., 3., 6., 8., 10.};
            break;
        case 1:
            non_ghost_ids  = {4, 5, 6, 7, 8};
            non_ghost_vals = {4., 5., 6., 7., 8.};
            ghost_ids      = {0, 9};
            expected       = {4., 5., 6., 7., 8., 0., 9.};
            break;
        case 2:
            non_ghost_ids  = {9, 10, 11};
            non_ghost_vals = {9., 10., 11.};
            ghost_ids      = {3, 5};
            expected       = {9., 10., 11., 3., 5.};
            break;
    }
    T_VECTOR x_with_ghosts(local_vec_size, ghost_ids, is_global_size);
    x_with_ghosts.set(non_ghost_ids, non_ghost_vals);
    MathLib::finalizeVectorAssembly(x_with_ghosts);

    ASSERT_EQ(12u, x_with_ghosts.size());

    std::vector<double> loc_v1(  x_with_ghosts.getLocalSize()
                               + x_with_ghosts.getGhostSize() );
    x_with_ghosts.copyValues(loc_v1);
    for (std::size_t i=0; i<expected.size(); i++)
    {
         ASSERT_EQ(expected[i], loc_v1[i]);
    }
}
#endif

} // end namespace

//--------------------------------------------
#if defined(USE_PETSC)
TEST(MPITest_Math, CheckInterface_PETScVector)
{
    checkGlobalVectorInterfacePETSc<MathLib::PETScVector >();
}
#elif defined(OGS_USE_EIGEN)
TEST(Math, CheckInterface_EigenVector)
{
    checkGlobalVectorInterface<MathLib::EigenVector >();
}
#endif
