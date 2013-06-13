/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Implementation tests.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>
#include <boost/property_tree/ptree.hpp>

#include "MathLib/LinAlg/Dense/DenseVector.h"
#include "MathLib/LinAlg/Dense/DenseMatrix.h"
#include "MathLib/LinAlg/Dense/GlobalDenseMatrix.h"
#include "MathLib/LinAlg/Dense/DenseTools.h"
#include "MathLib/LinAlg/FinishMatrixAssembly.h"
#include "MathLib/LinAlg/Solvers/DenseDirectLinearSolver.h"
#ifdef USE_LIS
#include "MathLib/LinAlg/Lis/LisVector.h"
#include "MathLib/LinAlg/Lis/LisMatrix.h"
#include "MathLib/LinAlg/Lis/LisLinearSolver.h"
#include "MathLib/LinAlg/Lis/LisTools.h"
#endif
#include "MathLib/LinAlg/Sparse/Sparsity.h"

#include "../TestTools.h"

namespace
{

template<class T_Mat>
void setMatrix9x9(T_Mat &mat)
{
    double d_mat[] = {
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
    for (unsigned i=0; i<9; i++)
        for (unsigned j=0; j<9; j++)
            if (d_mat[i*9+j]!=.0)
                mat.setValue(i, j, d_mat[i*9+j]);

}

struct Example1
{
    MathLib::GlobalDenseMatrix<double> mat;
    std::vector<size_t> list_dirichlet_bc_id;
    std::vector<double> list_dirichlet_bc_value;
    static const size_t dim_eqs = 9;
    std::vector<double> exH;
    MathLib::RowMajorSparsity sparse;

    Example1()
    : mat(dim_eqs, dim_eqs)
    {
        setMatrix9x9(mat);
        size_t int_dirichlet_bc_id[] = {2,5,8,0,3,6};
        list_dirichlet_bc_id.assign(int_dirichlet_bc_id, int_dirichlet_bc_id+6);
        list_dirichlet_bc_value.resize(6);
        fill(list_dirichlet_bc_value.begin(), list_dirichlet_bc_value.begin()+3, .0);
        fill(list_dirichlet_bc_value.begin()+3, list_dirichlet_bc_value.end(), 1.0);
        exH.resize(9);
        for (size_t i=0; i<9; i++) {
            if (i%3==0) exH[i] = 1.0;
            if (i%3==1) exH[i] = 0.5;
            if (i%3==2) exH[i] = 0.;
        }
        sparse.resize(dim_eqs);
        for (size_t i=0; i<dim_eqs; i++) {
            for (size_t j=0; j<dim_eqs; j++) {
                if (mat(i,j)!=.0)
                    sparse[i].insert(j);
            }
        }

    }
};

} // end namespace

TEST(Math, LinAlgDenseLinearSolver1)
{
    MathLib::GlobalDenseMatrix<double> A(9, 9);
    setMatrix9x9(A);
    finishMatrixAssembly(A);

    MathLib::DenseVector<double> b(A.getNRows());
    MathLib::DenseVector<double> x(A.getNRows());
    x = 1.0;
    A.matvec(x, b);
    x = 0.0;
    MathLib::DenseDirectLinearSolver ls;
    ls.solve(A, b, x);

    std::vector<double> exp_x(9, 1.0);
    ASSERT_DOUBLE_ARRAY_EQ(exp_x, x, 9);
}

TEST(Math, LinAlgDenseLinearSolver2)
{
    Example1 ex1;

    MathLib::GlobalDenseMatrix<double> A(ex1.dim_eqs, ex1.dim_eqs);
    MathLib::DenseVector<double> b(A.getNRows());
    MathLib::DenseVector<double> x(A.getNRows());

    for (size_t i=0; i<ex1.dim_eqs; i++)
        for (size_t j=0; j<ex1.dim_eqs; j++)
            A.addValue(i, j, ex1.mat(i, j));

    // apply BC
    MathLib::applyKnownSolution(A, b, ex1.list_dirichlet_bc_id, ex1.list_dirichlet_bc_value);
    finishMatrixAssembly(A);

    // solve
    MathLib::DenseDirectLinearSolver ls;
    ls.solve(A, b, x);

    ASSERT_DOUBLE_ARRAY_EQ(ex1.exH, x, ex1.dim_eqs, 1.e-5);
}


#ifdef USE_LIS
TEST(Math, LinAlgLisMatrix)
{
    MathLib::LisMatrix mat(9);
    setMatrix9x9(mat);
    ASSERT_EQ(9u, mat.getNRows());
    ASSERT_EQ(2.66667e-011, mat.getMaxDiagCoeff());
}

TEST(Math, LinAlgLisLinearSolver1)
{
    MathLib::LisMatrix A(9);
    setMatrix9x9(A);
    finishMatrixAssembly(A);

    MathLib::LisVector b(A.getNRows());
    MathLib::LisVector x(A.getNRows());
    x = 1.0;
    A.matvec(x, b);
    x = 0.0;
    MathLib::LisLinearSolver ls;
    ls.solve(A, b, x);

    std::vector<double> exp_x(9, 1.0);
    ASSERT_DOUBLE_ARRAY_EQ(exp_x, x, 9);
}

TEST(Math, LinAlgLisLinearSolver2)
{
    Example1 ex1;

    MathLib::LisMatrix A(ex1.dim_eqs);
    MathLib::LisVector b(A.getNRows());
    MathLib::LisVector x(A.getNRows());

    for (size_t i=0; i<ex1.dim_eqs; i++) {
        for (size_t j=0; j<ex1.dim_eqs; j++) {
            double v = ex1.mat(i, j);
            if (v!=.0)
                A.addValue(i, j, v);
        }
    }

    // apply BC
    MathLib::applyKnownSolution(A, b, ex1.list_dirichlet_bc_id, ex1.list_dirichlet_bc_value);
    finishMatrixAssembly(A);

    // set solver options using Boost property tree
    boost::property_tree::ptree t_root;
    boost::property_tree::ptree t_solver;
    t_solver.put("solver_type", "CG");
    t_solver.put("precon_type", "NONE");
    t_solver.put("matrix_type", "CCS");
    t_solver.put("error_tolerance", 1e-15);
    t_solver.put("max_iteration_step", 1000);
    t_root.put_child("LinearSolver", t_solver);
    MathLib::LisLinearSolver ls;
    ls.setOption(t_root);

    // check if the option was correctly parsed
    MathLib::LisOption &lisOption = ls.getOption();
    ASSERT_EQ(MathLib::LisOption::SolverType::CG, lisOption.solver_type);
    ASSERT_EQ(MathLib::LisOption::PreconType::NONE, lisOption.precon_type);
    ASSERT_EQ(MathLib::LisOption::MatrixType::CCS, lisOption.matrix_type);
    ASSERT_EQ(1e-15, lisOption.error_tolerance);
    ASSERT_EQ(1000, lisOption.max_iterations);

    // solve
    ls.solve(A, b, x);

    ASSERT_DOUBLE_ARRAY_EQ(ex1.exH, x, ex1.dim_eqs, 1.e-5);
}
#endif


