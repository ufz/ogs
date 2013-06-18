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

#include "MathLib/LinAlg/Dense/DenseVector.h"
#include "MathLib/LinAlg/Dense/DenseMatrix.h"
#include "MathLib/LinAlg/Dense/GlobalDenseMatrix.h"
#include "MathLib/LinAlg/Dense/DenseTools.h"
#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"

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

} // end namespace

TEST(MathLib, ApplyKnownSolutionAndSolveLinearSystem)
{
    Example1 ex1;

    MathLib::DenseVector<double> rhs(ex1.mat.getNRows());
    MathLib::DenseVector<double> x(ex1.mat.getNRows());

    // apply BC
    MathLib::applyKnownSolution(ex1.mat, rhs, ex1.vec_dirichlet_bc_id, ex1.vec_dirichlet_bc_value);

    // solve
    MathLib::GaussAlgorithm<MathLib::GlobalDenseMatrix<double>, MathLib::DenseVector<double> > ls(ex1.mat);
    ls.execute(rhs, x);

    ASSERT_DOUBLE_ARRAY_EQ(ex1.exH, x, ex1.mat.getNRows(), 1e-5);
}


