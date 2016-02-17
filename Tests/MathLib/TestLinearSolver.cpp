/**
 * \file
 * \author Norihiro Watanabe
 * \author Wenqing Wang
 * \date   2013-04-16, 2014-04
 * \brief  Implementation tests.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include "MathLib/LinAlg/Dense/DenseVector.h"
#include "MathLib/LinAlg/Dense/DenseMatrix.h"
#include "MathLib/LinAlg/Dense/GlobalDenseMatrix.h"
#include "MathLib/LinAlg/FinalizeMatrixAssembly.h"
#include "MathLib/LinAlg/ApplyKnownSolution.h"
#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"

#ifdef OGS_USE_EIGEN
#include "MathLib/LinAlg/Eigen/EigenMatrix.h"
#include "MathLib/LinAlg/Eigen/EigenVector.h"
#include "MathLib/LinAlg/Eigen/EigenLinearSolver.h"
#endif

#if defined(OGS_USE_EIGEN) && defined(USE_LIS)
#include "MathLib/LinAlg/EigenLis/EigenLisLinearSolver.h"
#endif

#ifdef USE_LIS
#include "MathLib/LinAlg/Lis/LisLinearSolver.h"
#endif

#ifdef USE_PETSC
#include "MathLib/LinAlg/PETSc/PETScMatrix.h"
#include "MathLib/LinAlg/PETSc/PETScVector.h"
#include "MathLib/LinAlg/PETSc/PETScLinearSolver.h"
#endif

#include "../TestTools.h"

namespace
{

template<class T_Mat>
void setMatrix9x9(T_Mat &mat)
{
    double d_mat[] =
    {
        6.66667e-012, -1.66667e-012, 0, -1.66667e-012, -3.33333e-012, 0, 0, 0, 0,
        -1.66667e-012, 1.33333e-011, -1.66667e-012, -3.33333e-012, -3.33333e-012, -3.33333e-012, 0, 0, 0,
        0, -1.66667e-012, 6.66667e-012, 0, -3.33333e-012, -1.66667e-012, 0, 0, 0,
        -1.66667e-012, -3.33333e-012, 0, 1.33333e-011, -3.33333e-012, 0, -1.66667e-012, -3.33333e-012, 0,
        -3.33333e-012, -3.33333e-012, -3.33333e-012, -3.33333e-012, 2.66667e-011, -3.33333e-012, -3.33333e-012, -3.33333e-012, -3.33333e-012,
        0, -3.33333e-012, -1.66667e-012, 0, -3.33333e-012, 1.33333e-011, 0, -3.33333e-012, -1.66667e-012,
        0, 0, 0, -1.66667e-012, -3.33333e-012, 0, 6.66667e-012, -1.66667e-012, 0,
        0, 0, 0, -3.33333e-012, -3.33333e-012, -3.33333e-012, -1.66667e-012, 1.33333e-011, -1.66667e-012,
        0, 0, 0, 0, -3.33333e-012, -1.66667e-012, 0, -1.66667e-012, 6.66667e-012
    };
    for (unsigned i = 0; i < 9; i++)
        for (unsigned j = 0; j < 9; j++)
            mat.setValue(i, j, d_mat[i*9+j]);
}

template<typename IntType> struct Example1
{
    MathLib::GlobalDenseMatrix<double> mat;
    std::vector<IntType> vec_dirichlet_bc_id;
    std::vector<double> vec_dirichlet_bc_value;
    static const std::size_t dim_eqs = 9;
    double* exH;

    Example1()
        : mat(dim_eqs, dim_eqs), exH(new double[dim_eqs])
    {
        setMatrix9x9(mat);
        IntType int_dirichlet_bc_id[] = {2,5,8,0,3,6};
        vec_dirichlet_bc_id.assign(int_dirichlet_bc_id, int_dirichlet_bc_id+6);
        vec_dirichlet_bc_value.resize(6);
        std::fill(vec_dirichlet_bc_value.begin(), vec_dirichlet_bc_value.begin()+3, .0);
        std::fill(vec_dirichlet_bc_value.begin()+3, vec_dirichlet_bc_value.end(), 1.0);
        for (std::size_t i=0; i<9; i++)
        {
            if (i%3==0) exH[i] = 1.0;
            if (i%3==1) exH[i] = 0.5;
            if (i%3==2) exH[i] = 0.;
        }
    }

    ~Example1()
    {
        delete [] exH;
    }
};

template <class T_MATRIX, class T_VECTOR, class T_LINEAR_SOVLER, typename IntType>
void checkLinearSolverInterface(T_MATRIX &A, BaseLib::ConfigTree const& ls_option)
{
    Example1<IntType> ex1;

    // set a coefficient matrix
    A.setZero();
    for (std::size_t i=0; i<ex1.dim_eqs; i++)
    {
        for (std::size_t j=0; j<ex1.dim_eqs; j++)
        {
            double v = ex1.mat(i, j);
            if (v!=.0)
                A.add(i, j, v);
        }
    }

    // set RHS and solution vectors
    T_VECTOR rhs(ex1.dim_eqs);
    T_VECTOR x(ex1.dim_eqs);

    // apply BC
    MathLib::applyKnownSolution(A, rhs, x, ex1.vec_dirichlet_bc_id, ex1.vec_dirichlet_bc_value);

    MathLib::finalizeMatrixAssembly(A);

    // solve
    T_LINEAR_SOVLER ls(A, "dummy_name", &ls_option);
    ls.solve(rhs, x);

    ASSERT_ARRAY_NEAR(ex1.exH, x, ex1.dim_eqs, 1e-5);

}

#ifdef USE_PETSC
template <class T_MATRIX, class T_VECTOR, class T_LINEAR_SOVLER>
void checkLinearSolverInterface(T_MATRIX& A, T_VECTOR& b,
                                const std::string& prefix_name,
                                BaseLib::ConfigTree const& ls_option)
{
    int mrank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &mrank);
    // Add entries
    MathLib::DenseMatrix<double> loc_m(2, 2);
    loc_m(0, 0) = 1. +  mrank;
    loc_m(0, 1) = 2. +  mrank;
    loc_m(1, 0) = 3. +  mrank;
    loc_m(1, 1) = 4. +  mrank;

    std::vector<int> row_pos(2);
    std::vector<int> col_pos(2);
    row_pos[0] = 2 * mrank;
    row_pos[1] = 2 * mrank + 1;
    col_pos[0] = row_pos[0];
    col_pos[1] = row_pos[1];

    A.add(row_pos, col_pos, loc_m);

    MathLib::finalizeMatrixAssembly(A);

    const bool deep_copy = false;
    T_VECTOR x(b, deep_copy);

    std::vector<double> local_vec(2);
    local_vec[0] = mrank+1;
    local_vec[1] = 2. * (mrank+1);
    x.set(row_pos, local_vec);

    double x0[6];
    double x1[6];
    x.getGlobalVector(x0);

    A.multiply(x, b);

    // apply BC
    std::vector<int> bc_id;  // Type must be int to match Petsc_Int
    std::vector<double> bc_value;

    if(mrank == 1)
    {
        bc_id.resize(1);
        bc_value.resize(1);
        bc_id[0] = 2 * mrank;
        bc_value[0] = mrank+1;
    }

    MathLib::applyKnownSolution(A, b, x, bc_id, bc_value);

    MathLib::finalizeMatrixAssembly(A);

    // solve
    T_LINEAR_SOVLER ls(A, prefix_name, &ls_option);
    EXPECT_TRUE(ls.solve(b, x));

    EXPECT_GT(ls.getNumberOfIterations(), 0u);

    x.getGlobalVector(x1);
    ASSERT_ARRAY_NEAR(x0, x1, 6, 1e-5);
}
#endif

} // end namespace

TEST(MathLib, CheckInterface_GaussAlgorithm)
{
    boost::property_tree::ptree t_root;
    BaseLib::ConfigTree conf(t_root, "",
        BaseLib::ConfigTree::onerror, BaseLib::ConfigTree::onwarning);

    using Example = Example1<std::size_t>;

    typedef MathLib::GaussAlgorithm<MathLib::GlobalDenseMatrix<double>, MathLib::DenseVector<double> > LinearSolverType;
    MathLib::GlobalDenseMatrix<double> A(Example::dim_eqs, Example::dim_eqs);
    checkLinearSolverInterface<MathLib::GlobalDenseMatrix<double>,
                               MathLib::DenseVector<double>, LinearSolverType, std::size_t>(
        A, conf);
}

#ifdef OGS_USE_EIGEN
TEST(Math, CheckInterface_Eigen)
{
    // set solver options using Boost property tree
    boost::property_tree::ptree t_root;
    boost::property_tree::ptree t_solver;
    t_solver.put("solver_type", "CG");
    t_solver.put("precon_type", "NONE");
    t_solver.put("error_tolerance", 1e-15);
    t_solver.put("max_iteration_step", 1000);
    t_root.put_child("eigen", t_solver);
    BaseLib::ConfigTree conf(t_root, "",
        BaseLib::ConfigTree::onerror, BaseLib::ConfigTree::onwarning);

    using IntType = MathLib::EigenMatrix::IndexType;

    MathLib::EigenMatrix A(Example1<IntType>::dim_eqs);
    checkLinearSolverInterface<MathLib::EigenMatrix, MathLib::EigenVector,
                               MathLib::EigenLinearSolver, IntType>(A, conf);
}
#endif

#if defined(OGS_USE_EIGEN) && defined(USE_LIS)
TEST(Math, CheckInterface_EigenLis)
{
    // set solver options using Boost property tree
    boost::property_tree::ptree t_root;
    boost::property_tree::ptree t_solver;
    t_root.put("lis", "-i cg -p none -tol 1e-15 -maxiter 1000");
    BaseLib::ConfigTree conf(t_root, "",
        BaseLib::ConfigTree::onerror, BaseLib::ConfigTree::onwarning);

    using IntType = MathLib::LisMatrix::IndexType;

    MathLib::EigenMatrix A(Example1<IntType>::dim_eqs);
    checkLinearSolverInterface<MathLib::EigenMatrix, MathLib::EigenVector,
                               MathLib::EigenLisLinearSolver, IntType>(A, conf);
}
#endif

#ifdef USE_LIS
TEST(Math, CheckInterface_Lis)
{
    // set solver options using Boost property tree
    boost::property_tree::ptree t_root;
    boost::property_tree::ptree t_solver;
    t_root.put("lis", "-i cg -p none -tol 1e-15 -maxiter 1000");
    BaseLib::ConfigTree conf(t_root, "",
        BaseLib::ConfigTree::onerror, BaseLib::ConfigTree::onwarning);

    using IntType = MathLib::LisMatrix::IndexType;

    MathLib::LisMatrix A(Example1<IntType>::dim_eqs);
    checkLinearSolverInterface<MathLib::LisMatrix, MathLib::LisVector,
                               MathLib::LisLinearSolver, IntType>(A, conf);
}
#endif

#ifdef USE_PETSC
TEST(MPITest_Math, CheckInterface_PETSc_Linear_Solver_basic)
{
    MathLib::PETScMatrixOption opt;
    opt.d_nz = 2;
    opt.o_nz = 0;
    opt.is_global_size = false;
    opt.n_local_cols = 2;
    MathLib::PETScMatrix A(2, opt);

    const bool is_global_size = false;
    MathLib::PETScVector b(2, is_global_size);

    boost::property_tree::ptree t_root;
    t_root.put("petsc",
                 "-ptest1_ksp_type bcgs "
                 "-ptest1_ksp_rtol 1.e-8 "
                 "-ptest1_ksp_atol 1.e-50 "
                 "-ptest1_ksp_max_it 1000 "
                 "-ptest1_pc_type bjacobi");

    checkLinearSolverInterface<MathLib::PETScMatrix,
                               MathLib::PETScVector,
                               MathLib::PETScLinearSolver>(
        A, b, "ptest1_",
        BaseLib::ConfigTree(t_root, "", BaseLib::ConfigTree::onerror,
                            BaseLib::ConfigTree::onwarning)
    );
}

TEST(MPITest_Math, CheckInterface_PETSc_Linear_Solver_chebyshev_sor)
{
    MathLib::PETScMatrixOption opt;
    opt.d_nz = 2;
    opt.o_nz = 0;
    opt.is_global_size = false;
    opt.n_local_cols = 2;
    MathLib::PETScMatrix A(2, opt);

    const bool is_global_size = false;
    MathLib::PETScVector b(2, is_global_size);

    boost::property_tree::ptree t_root;
    t_root.put("petsc",
                 "-ptest2_ksp_type chebyshev "
                 "-ptest2_ksp_rtol 1.e-8 "
                 "-ptest2_ksp_atol 1.e-50"
                 "-ptest2_ksp_max_it 1000 "
                 "-ptest2_pc_type sor");

    checkLinearSolverInterface<MathLib::PETScMatrix,
                               MathLib::PETScVector,
                               MathLib::PETScLinearSolver>(
        A, b, "ptest2_",
        BaseLib::ConfigTree(t_root, "", BaseLib::ConfigTree::onerror,
                            BaseLib::ConfigTree::onwarning)
    );
}

TEST(MPITest_Math, CheckInterface_PETSc_Linear_Solver_gmres_amg)
{
    MathLib::PETScMatrixOption opt;
    opt.d_nz = 2;
    opt.o_nz = 0;
    opt.is_global_size = false;
    opt.n_local_cols = 2;
    MathLib::PETScMatrix A(2, opt);

    const bool is_global_size = false;
    MathLib::PETScVector b(2, is_global_size);

    boost::property_tree::ptree t_root;
    t_root.put("petsc",
                 "-ptest3_ksp_type gmres "
                 "-ptest3_ksp_rtol 1.e-8 "
                 "-ptest3_ksp_gmres_restart 20 "
                 "-ptest3_ksp_gmres_classicalgramschmidt "
                 "-ptest3_pc_type gamg "
                 "-ptest3_pc_gamg_type agg "
                 "-ptest3_pc_gamg_agg_nsmooths 2");

    checkLinearSolverInterface<MathLib::PETScMatrix,
                               MathLib::PETScVector,
                               MathLib::PETScLinearSolver>(
        A, b, "ptest3_",
        BaseLib::ConfigTree(t_root, "", BaseLib::ConfigTree::onerror,
                            BaseLib::ConfigTree::onwarning)
    );
}

#endif


