/**
 * @file TestDenseGaussAlgorithm.cpp
 * @author Thomas Fischer
 * @date Jun 17, 2013
 * @brief
 *
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <cstdlib>
#include <ctime>
#include <limits>
#include <algorithm>

#include <gtest/gtest.h>

#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"

TEST(MathLib, DenseGaussAlgorithm)
{
    std::size_t n_rows(100);
    std::size_t n_cols(n_rows);

    MathLib::DenseMatrix<double,std::size_t> mat(n_rows, n_cols);

    // *** fill matrix with arbitrary values
    // ** initialize random seed
    srand ( static_cast<unsigned>(time(nullptr)) );
    // ** loop over rows and columns
    for (std::size_t i(0); i<n_rows; i++) {
        for (std::size_t j(0); j<n_cols; j++) {
            mat(i,j) = rand()/static_cast<double>(RAND_MAX);
        }
    }

    // *** create solution vector, set all entries to 0.0
    double *x(new double[n_cols]);
    std::fill(x,x+n_cols, 0.0);
    double *b0(mat * x);

    // *** create other right hand sides,
    // set all entries of the solution vector to 1.0
    std::fill(x,x+n_cols, 1.0);
    double *b1(mat * x);

    std::generate(x,x+n_cols, std::rand);
    double *b2(mat * x);

    // right hand side and solution vector with random entries
    double *b3(mat * x);
    double *b3_copy(mat * x);
    double *x3 (new double[n_cols]);
    std::generate(x3,x3+n_cols, std::rand);

    MathLib::GaussAlgorithm<MathLib::DenseMatrix<double, std::size_t>, double*> gauss;

    // solve with b0 as right hand side
    gauss.solve(mat, b0, true);
    for (std::size_t i(0); i<n_rows; i++) {
        ASSERT_NEAR(b0[i], 0.0, std::numeric_limits<float>::epsilon());
    }

    // solve with b1 as right hand side
    gauss.solve(mat, b1, false);
    for (std::size_t i(0); i<n_rows; i++) {
        ASSERT_NEAR(b1[i], 1.0, std::numeric_limits<float>::epsilon());
    }

    // solve with b2 as right hand side
    gauss.solve(mat, b2, false);
    for (std::size_t i(0); i<n_rows; i++) {
        ASSERT_NEAR(fabs(b2[i]-x[i])/fabs(x[i]), 0.0, std::numeric_limits<float>::epsilon());
    }

    // solve with b3 as right hand side and x3 as solution vector
    gauss.solve(mat, b3, x3, false);
    for (std::size_t i(0); i<n_rows; i++) {
        ASSERT_NEAR(fabs(x3[i]-x[i])/fabs(x[i]), 0.0, std::numeric_limits<float>::epsilon());
    }
    // assure entries of vector b3 are not changed
    for (std::size_t i(0); i<n_rows; i++) {
        ASSERT_NEAR(fabs(b3[i]-b3_copy[i])/fabs(b3[i]), 0.0, std::numeric_limits<float>::epsilon());
    }

    delete [] x;
    delete [] b0;
    delete [] b1;
    delete [] b2;
    delete [] b3;
    delete [] x3;
    delete [] b3_copy;
}

