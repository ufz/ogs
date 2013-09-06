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
#include "MathLib/LinAlg/FinalizeMatrixAssembly.h"
#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"
#ifdef USE_LIS
#include "MathLib/LinAlg/Lis/LisLinearSolver.h"
#include "MathLib/LinAlg/Lis/LisTools.h"
#endif

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
	for (unsigned i = 0; i < 9; i++)
		for (unsigned j = 0; j < 9; j++)
			mat.setValue(i, j, d_mat[i*9+j]);
}

struct Example1
{
    MathLib::GlobalDenseMatrix<double> mat;
    std::vector<size_t> vec_dirichlet_bc_id;
    std::vector<double> vec_dirichlet_bc_value;
    static const std::size_t dim_eqs = 9;
    double* exH;

    Example1()
    : mat(dim_eqs, dim_eqs), exH(new double[dim_eqs])
    {
        setMatrix9x9(mat);
        std::size_t int_dirichlet_bc_id[] = {2,5,8,0,3,6};
        vec_dirichlet_bc_id.assign(int_dirichlet_bc_id, int_dirichlet_bc_id+6);
        vec_dirichlet_bc_value.resize(6);
        std::fill(vec_dirichlet_bc_value.begin(), vec_dirichlet_bc_value.begin()+3, .0);
        std::fill(vec_dirichlet_bc_value.begin()+3, vec_dirichlet_bc_value.end(), 1.0);
        for (std::size_t i=0; i<9; i++) {
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

template <class T_MATRIX, class T_VECTOR, class T_LINEAR_SOVLER>
void checkLinearSolverInterface(T_MATRIX &A, boost::property_tree::ptree &ls_option)
{
    Example1 ex1;

    // set a coefficient matrix
    A.setZero();
    for (size_t i=0; i<ex1.dim_eqs; i++) {
        for (size_t j=0; j<ex1.dim_eqs; j++) {
            double v = ex1.mat(i, j);
            if (v!=.0)
                A.add(i, j, v);
        }
    }

    // set RHS and solution vectors
    T_VECTOR rhs(ex1.dim_eqs);
    T_VECTOR x(ex1.dim_eqs);

    // apply BC
    MathLib::applyKnownSolution(A, rhs, ex1.vec_dirichlet_bc_id, ex1.vec_dirichlet_bc_value);

    MathLib::finalizeMatrixAssembly(A);

    // solve
    T_LINEAR_SOVLER ls(A, &ls_option);
    ls.solve(rhs, x);

    ASSERT_ARRAY_NEAR(ex1.exH, x, ex1.dim_eqs, 1e-5);

}

} // end namespace

TEST(MathLib, CheckInterface_GaussAlgorithm)
{
    boost::property_tree::ptree t_root;
    boost::property_tree::ptree t_solver;
    t_root.put_child("LinearSolver", t_solver);

    typedef MathLib::GaussAlgorithm<MathLib::GlobalDenseMatrix<double>, MathLib::DenseVector<double> > LinearSolverType;
    MathLib::GlobalDenseMatrix<double> A(Example1::dim_eqs, Example1::dim_eqs);
    checkLinearSolverInterface<MathLib::GlobalDenseMatrix<double>, MathLib::DenseVector<double>, LinearSolverType>(A, t_root);
}

#ifdef USE_LIS
TEST(Math, CheckInterface_Lis)
{
    // set solver options using Boost property tree
    boost::property_tree::ptree t_root;
    boost::property_tree::ptree t_solver;
    t_solver.put("solver_type", "CG");
    t_solver.put("precon_type", "NONE");
    t_solver.put("error_tolerance", 1e-15);
    t_solver.put("max_iteration_step", 1000);
    t_root.put_child("LinearSolver", t_solver);

    MathLib::LisMatrix A(Example1::dim_eqs);
    checkLinearSolverInterface<MathLib::LisMatrix, MathLib::LisVector, MathLib::LisLinearSolver>(A, t_root);
}
#endif

