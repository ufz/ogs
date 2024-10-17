/**
 * \file
 * \author Norihiro Watanabe
 * \author Wenqing Wang
 * \date   2013-05-15 -- 2014-01
 * \brief  Interface tests of global vector classes
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include <cmath>
#include <numeric>
#include <ranges>

#include "MathLib/LinAlg/LinAlg.h"
#include "Tests/TestTools.h"

#if defined(USE_PETSC)
#include "MathLib/LinAlg/PETSc/PETScVector.h"
#else
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
    x.setZero();

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
    // y -= x;
    axpy(y, -1., x);
    ASSERT_EQ(2.0, y.get(0));
    // y = 1.0;
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
double getEuclideanNorm(std::vector<double> const& v)
{
    // Create a view that squares each element
    auto squared_view =
        v | std::views::transform([](double x) { return x * x; });

    return std::sqrt(
        std::accumulate(squared_view.begin(), squared_view.end(), 0.0));
}

template <class T_VECTOR>
void checkPETScVectorNoExplictGhostID(T_VECTOR& x, int const msize)
{
    ASSERT_EQ(x.getRangeEnd() - x.getRangeBegin(), x.getLocalSize());

    // Note id0_rank may vary because of the PETSc determined partitioning
    const int id0_rank = x.getRangeBegin();
    x.setLocalAccessibleVector();
    ASSERT_EQ(.0, x.get(id0_rank));

    set(x, 10.);  // x:=10.0
    x.finalizeAssembly();

    // Value of x is not copied to y
    const bool deep_copy = false;
    T_VECTOR y(x, deep_copy);  // y:= 0.0
    y.setLocalAccessibleVector();
    ASSERT_EQ(0, y.get(id0_rank));

    set(y, 10.0);  // y:=10.0
    y.setLocalAccessibleVector();
    ASSERT_EQ(10, y.get(id0_rank));

    // y += x
    axpy(y, 1., x);  // y = x + y := 20
    y.setLocalAccessibleVector();
    ASSERT_EQ(20, y.get(id0_rank));

    // y -= x
    axpy(y, -1., x);  // y = y - x : = 10.0
    y.setLocalAccessibleVector();
    ASSERT_EQ(10, y.get(id0_rank));

    // Set values by each rank
    std::vector<double> local_vec(2, 10.0);
    std::vector<GlobalIndexType> vec_pos(2);

    vec_pos[0] = id0_rank;
    int const id0_pls_1 = id0_rank + 1;
    // If id0_pls_1 is out of the current rank, set its negative value to
    // vec_pos[1] to avoid adding value.
    vec_pos[1] = id0_pls_1 < y.getRangeEnd() ? id0_pls_1 : -id0_pls_1;
    y.add(vec_pos, local_vec);
    y.finalizeAssembly();

    int const x_size = x.size();
    // Count the number of the change elements of y by y.add
    std::vector<int> all_vec_pos(msize * 2);
    MPI_Allgather(vec_pos.data(), 2, MPI_INT, all_vec_pos.data(), 2, MPI_INT,
                  PETSC_COMM_WORLD);

    std::vector<double> x_raw(x_size);
    x.getGlobalVector(x_raw);
    std::ranges::for_each(
        all_vec_pos,
        [&x_raw](int id)
        {
            if (id >= 0 && id < static_cast<int>(x_raw.size()))
            {
                x_raw[id] = 20.0;
            }
        });

    std::vector<double> y_raw(x_size);

    y.getGlobalVector(y_raw);

    ASSERT_ARRAY_NEAR(x_raw, y_raw, x_size, 1e-10);

    ASSERT_NEAR(getEuclideanNorm(y_raw), norm2(y), 1e-10);

    // Deep copy
    T_VECTOR x_deep_copied(y);
    x_deep_copied.getGlobalVector(y_raw);
    ASSERT_ARRAY_NEAR(x_raw, y_raw, x_size, 1e-10);
}

template <class T_VECTOR>
void checkPETScVectorExplictGhostID()
{
    int msize;
    MPI_Comm_size(PETSC_COMM_WORLD, &msize);

    ASSERT_EQ(3u, msize);

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
    int mrank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &mrank);
    std::size_t local_vec_size = 4;
    if (mrank == 1)
    {
        local_vec_size = 5;
    }
    else if (mrank == 2)
    {
        local_vec_size = 3;
    }
    std::vector<GlobalIndexType> non_ghost_ids(local_vec_size);
    std::vector<double> non_ghost_vals(local_vec_size);
    std::size_t nghosts = 3;
    if (mrank > 0)
    {
        nghosts = 2;
    }
    std::vector<GlobalIndexType> ghost_ids(nghosts);
    std::vector<double> expected;
    switch (mrank)
    {
        case 0:
            non_ghost_ids = {0, 1, 2, 3};
            non_ghost_vals = {0., 1., 2., 3.};
            ghost_ids = {6, 8, 10};
            expected = {0., 1., 2., 3., 6., 8., 10.};
            break;
        case 1:
            non_ghost_ids = {4, 5, 6, 7, 8};
            non_ghost_vals = {4., 5., 6., 7., 8.};
            ghost_ids = {0, 9};
            expected = {4., 5., 6., 7., 8., 0., 9.};
            break;
        case 2:
            non_ghost_ids = {9, 10, 11};
            non_ghost_vals = {9., 10., 11.};
            ghost_ids = {3, 5};
            expected = {9., 10., 11., 3., 5.};
            break;
    }
    T_VECTOR x_with_ghosts(local_vec_size, ghost_ids, false /*is_global_size*/);
    x_with_ghosts.set(non_ghost_ids, non_ghost_vals);
    MathLib::finalizeVectorAssembly(x_with_ghosts);

    ASSERT_EQ(12u, x_with_ghosts.size());

    std::vector<double> loc_v1(x_with_ghosts.getLocalSize() +
                               x_with_ghosts.getGhostSize());
    x_with_ghosts.copyValues(loc_v1);
    for (std::size_t i = 0; i < expected.size(); i++)
    {
        ASSERT_EQ(expected[i], loc_v1[i]);
    }
}
#endif

}  // end namespace

//--------------------------------------------
#if defined(USE_PETSC)
TEST(MPITest_Math, PETScVectorPatitionedAutomatically)
{
    int msize;
    MPI_Comm_size(PETSC_COMM_WORLD, &msize);

    MathLib::PETScVector x(16);  // x:=0.0
    ASSERT_LE(msize, 16u);

    ASSERT_EQ(16u, x.size());
    checkPETScVectorNoExplictGhostID<MathLib::PETScVector>(x, msize);
}

TEST(MPITest_Math, PETScVectorFixedPartition)
{
    int msize;
    MPI_Comm_size(PETSC_COMM_WORLD, &msize);

    // Each partition has two elements.
    int const num_elements_per_rank = 2;
    MathLib::PETScVector x_fixed_p(num_elements_per_rank,
                                   false /*is_global_size*/);

    ASSERT_EQ(num_elements_per_rank * msize, x_fixed_p.size());

    int mrank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &mrank);
    ASSERT_EQ(2 * mrank, x_fixed_p.getRangeBegin());
    ASSERT_EQ(2 * mrank + 2, x_fixed_p.getRangeEnd());
    checkPETScVectorNoExplictGhostID<MathLib::PETScVector>(x_fixed_p, msize);
}

TEST(MPITest_Math, CheckPETScVectorExplictGhostID)
{
    checkPETScVectorExplictGhostID<MathLib::PETScVector>();
}
#else
TEST(Math, CheckInterface_EigenVector)
{
    checkGlobalVectorInterface<MathLib::EigenVector>();
}
#endif
