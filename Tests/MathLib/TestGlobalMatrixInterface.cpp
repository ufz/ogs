/**
 * \file
 * \author Norihiro Watanabe
 * \author Wenqing Wang
 * \date   2013-05-15, 2014-02
 * \brief  Interface tests of global matrix classes
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include "MathLib/LinAlg/Dense/DenseMatrix.h"
#include "MathLib/LinAlg/Dense/GlobalDenseMatrix.h"
#include "MathLib/LinAlg/FinalizeMatrixAssembly.h"

#ifdef USE_LIS
#include "MathLib/LinAlg/Lis/LisMatrix.h"
#endif

#ifdef USE_PETSC
#include "MathLib/LinAlg/PETSc/PETScMatrix.h"
#endif

namespace
{

template <class T_MATRIX>
void checkGlobalMatrixInterface(T_MATRIX &m)
{
    ASSERT_EQ(10u, m.getNRows());
    ASSERT_EQ(10u, m.getNCols());
    ASSERT_EQ(0u,  m.getRangeBegin());
    ASSERT_EQ(10u, m.getRangeEnd());

    m.setValue(0, 0, 1.0);
    m.add(0, 0, 1.0);
    m.setZero();

    MathLib::DenseMatrix<double> local_m(2, 2, 1.0);
    std::vector<std::size_t> vec_pos(2);
    vec_pos[0] = 1;
    vec_pos[1] = 3;
    m.add(vec_pos, vec_pos, local_m);

    ASSERT_TRUE(finalizeMatrixAssembly(m));
}

#ifdef USE_PETSC // or MPI
template <class T_MATRIX, class T_VECTOR>
void checkGlobalMatrixInterfaceMPI(T_MATRIX &m, T_VECTOR &v)
{
    int msize;
    MPI_Comm_size(PETSC_COMM_WORLD, &msize);
    int mrank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &mrank);

    ASSERT_EQ(3u, msize);
    ASSERT_EQ(m.getRangeEnd()-m.getRangeBegin(), m.getNLocalRows());

    int gathered_rows;
    int local_rows = m.getNLocalRows();
    MPI_Allreduce(&local_rows, &gathered_rows, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
    ASSERT_EQ(m.getNRows(), gathered_rows);

    int gathered_cols;
    int local_cols = m.getNLocalColumns();
    MPI_Allreduce(&local_cols, &gathered_cols, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
    ASSERT_EQ(m.getNCols(), gathered_cols);

    // Add entries
    MathLib::DenseMatrix<double> loc_m(2, 2);
    loc_m(0, 0) = 1.;
    loc_m(0, 1) = 2.;
    loc_m(1, 0) = 3.;
    loc_m(1, 1) = 4.;

    std::vector<int> row_pos(2);
    std::vector<int> col_pos(2);
    row_pos[0] = 2 * mrank;
    row_pos[1] = 2 * mrank + 1;
    col_pos[0] = row_pos[0];
    col_pos[1] = row_pos[1];

    m.add(row_pos, col_pos, loc_m);

    MathLib::finalizeMatrixAssembly(m);

    // Multiply by a vector
    v = 1.;
    const bool deep_copy = false;
    T_VECTOR y(v, deep_copy);
    m.multi(v, y);

    ASSERT_EQ(sqrt(3*(3*3 + 7*7)), y.getNorm());

    // Set/add a value
    if(mrank == 1)
    {
        m.set(0, 0, 5.0);
        m.add(1, 1, 5.0);
    }
    MathLib::finalizeMatrixAssembly(m);
    m.multi(v, y);

    ASSERT_EQ(sqrt((2*3*3 + 8*8 + 3*7*7)), y.getNorm());
}

// Rectanglular matrix
template <class T_MATRIX, class T_VECTOR>
void checkGlobalRectangularMatrixInterfaceMPI(T_MATRIX &m, T_VECTOR &v)
{
    int mrank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &mrank);

    ASSERT_EQ(m.getRangeEnd()-m.getRangeBegin(), m.getNLocalRows());

    int gathered_rows;
    int local_rows = m.getNLocalRows();
    MPI_Allreduce(&local_rows, &gathered_rows, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
    ASSERT_EQ(m.getNRows(), gathered_rows);

    int gathered_cols;
    int local_cols = m.getNLocalColumns();
    MPI_Allreduce(&local_cols, &gathered_cols, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
    ASSERT_EQ(m.getNCols(), gathered_cols);

    printf("m.getNRows()%d  m.getNRows() %d  %d",  m.getNRows(),  m.getNCols(), v.size() );

    // Add entries
    MathLib::DenseMatrix<double> loc_m(2, 3);
    loc_m(0, 0) = 1.;
    loc_m(0, 1) = 2.;
    loc_m(0, 1) = 3.;
    loc_m(1, 0) = 1.;
    loc_m(1, 1) = 2.;
    loc_m(1, 1) = 3.;

    std::vector<int> row_pos(2);
    std::vector<int> col_pos(3);
    row_pos[0] = 2 * mrank;
    row_pos[1] = 2 * mrank + 1;
    col_pos[0] = 3 * mrank;
    col_pos[1] = 3 * mrank + 1;
    col_pos[2] = 3 * mrank + 2;
    
    m.add(row_pos, col_pos, loc_m);

    MathLib::finalizeMatrixAssembly(m);

    /*
    Cannot perform the following multiplication due to PETSc bug.
    // Multiply by a vector
    v = 1.;
    const bool deep_copy = false;
    T_VECTOR y(v, deep_copy);
    m.multi(v, y);

    ASSERT_EQ(sqrt(24.), y.getNorm());
    */ 
}

#endif // end of: ifdef USE_PETSC // or MPI

} // end namespace

TEST(Math, CheckInterface_GlobalDenseMatrix)
{
    MathLib::GlobalDenseMatrix<double> m(10, 10);
    checkGlobalMatrixInterface(m);
}

#ifdef USE_LIS
TEST(Math, CheckInterface_LisMatrix)
{
    MathLib::LisMatrix m(10);
    checkGlobalMatrixInterface(m);
}
#endif

#ifdef USE_PETSC // or MPI
TEST(Math, CheckInterface_PETScMatrix_Local_Size)
{
    MathLib::PETScMatrixOption opt;
    opt._d_nz = 2;
    opt._o_nz = 0;
    opt._is_global_size = false;
    opt._n_local_cols = 2;
    MathLib::PETScMatrix A(2, opt);

    const bool is_gloabal_size = false;
    MathLib::PETScVector x(2, is_gloabal_size);

    checkGlobalMatrixInterfaceMPI(A, x);
}

TEST(Math, CheckInterface_PETScMatrix_Global_Size)
{
    MathLib::PETScMatrixOption opt;
    opt._d_nz = 2;
    opt._o_nz = 0;
    MathLib::PETScMatrix A(6, opt);

    MathLib::PETScVector x(6);

    checkGlobalMatrixInterfaceMPI(A, x);
}

// Test rectangular matrix
TEST(Math, CheckInterface_PETSc_Rectangular_Matrix_Local_Size)
{
    MathLib::PETScMatrixOption opt;
    opt._d_nz = 3;
    opt._o_nz = 0;
    opt._is_global_size = false;
    opt._n_local_cols = -1;
    MathLib::PETScMatrix A(2, 3, opt);

    const bool is_gloabal_size = false;
    MathLib::PETScVector x(3, is_gloabal_size);

    checkGlobalRectangularMatrixInterfaceMPI(A, x);
}

TEST(Math, CheckInterface_PETSc_Rectangular_Matrix_Global_Size)
{
    MathLib::PETScMatrixOption opt;
    opt._d_nz = 3;
    opt._o_nz = 0;
    MathLib::PETScMatrix A(6, 9, opt);

    MathLib::PETScVector x(9);

    checkGlobalRectangularMatrixInterfaceMPI(A, x);
}

#endif // end of: ifdef USE_PETSC // or MPI

