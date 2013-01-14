/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief  Implementation tests.
 *
 * \copyright
 * Copyright (c)  2013, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>
#include <boost/property_tree/ptree.hpp>

#include "MathLib/LinAlg/Sparse/Sparsity.h"
#ifdef USE_LIS
#include "MathLib/LinAlg/SystemOfLinearEquations/LisLinearSystem.h"
#endif

namespace
{

inline void ASSERT_DOUBLE_ARRAY_EQ(const double* Expected, const double* Actual, size_t N, double epsilon=1.0e-8) {
    for (size_t i=0; i<N; i++) \
        ASSERT_NEAR(Expected[i], Actual[i], epsilon);
}

struct Example1
{
    std::vector<double> mat;
    std::vector<size_t> list_dirichlet_bc_id;
    std::vector<double> list_dirichlet_bc_value;
    static const size_t dim_eqs = 9;
    std::vector<double> exH;
    MathLib::RowMajorSparsity sparse;

    Example1()
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
        mat.assign(d_mat, d_mat+dim_eqs*dim_eqs);
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
                if (mat[i*dim_eqs+j]!=.0)
                    sparse[i].insert(j);
            }
        }

    }
};
}

#ifdef USE_LIS
TEST(Math, LinearSystemLis)
{
    // set a problem
    Example1 ex1;

    // create a linear system
    MathLib::LisLinearSystem eqs(ex1.dim_eqs);

    // construct
    for (size_t i=0; i<ex1.dim_eqs; i++) {
        for (size_t j=0; j<ex1.dim_eqs; j++) {
            double v = ex1.mat[i*ex1.dim_eqs+j];
            if (v!=.0)
                eqs.addMatEntry(i, j, v);
        }
    }

    // apply BC
    eqs.setKnownSolution(ex1.list_dirichlet_bc_id, ex1.list_dirichlet_bc_value);

    // set solver options using Boost property tree
    boost::property_tree::ptree t_root;
    boost::property_tree::ptree t_solver;
    t_solver.put("solver_type", "CG");
    t_solver.put("precon_type", "NONE");
    t_solver.put("matrix_type", "CCS");
    t_solver.put("error_tolerance", 1e-15);
    t_solver.put("max_iteration_step", 1000);
    t_root.put_child("LinearSolver", t_solver);
    eqs.setOption(t_root);

    // check if the option was correctly parsed
    MathLib::LisOption &lisOption = eqs.getOption();
    ASSERT_EQ(MathLib::LisOption::SolverType::CG, lisOption.solver_type);
    ASSERT_EQ(MathLib::LisOption::PreconType::NONE, lisOption.precon_type);
    ASSERT_EQ(MathLib::LisOption::MatrixType::CCS, lisOption.matrix_type);
    ASSERT_EQ(1e-15, lisOption.error_tolerance);
    ASSERT_EQ(1000, lisOption.max_iterations);

    // solve
    eqs.solve();

//    eqs.printout();

    // check solution
    std::vector<double> vec_x(eqs.getDimension());
    eqs.getSolVec(&vec_x[0]);
    ASSERT_DOUBLE_ARRAY_EQ(&ex1.exH[0], &vec_x[0], ex1.dim_eqs, 1.e-5);
}

#endif

