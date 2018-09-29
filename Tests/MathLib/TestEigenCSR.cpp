/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include <vector>
#include <numeric>

#include <Eigen/Sparse>


/**
 * This test case checks if the internal Eigen::SparseMatrix compressed storage format
 * is a conventional CSR matrix. Currently this is the case, but it is not guaranteed
 * for all time.
 *
 * Cf. section "Sparse matrix format" on page
 * http://eigen.tuxfamily.org/dox/group__TutorialSparse.html
 */
TEST(MathLibEigen, Eigen2CSR)
{
    const int nrows = 16;
    const int ncols = nrows;

    Eigen::SparseMatrix<double, Eigen::RowMajor> mat(nrows, ncols);

    // set up sparsity pattern
    std::vector<int> pat(nrows);
    for (std::size_t i=0; i<nrows; ++i)
    {
        if (i==0 || i==nrows-1) {
            pat[i] = 2;
        } else {
            pat[i] = 3;
        }
    }

    // CSR representation of the matrix
    std::vector<double> values;
    std::vector<int> ia(nrows+1); // row offsets
    std::vector<int> ja;          // column indices

    std::partial_sum(pat.begin(), pat.end(), ia.begin()+1);

    const int nnz = ia.back();
    values.reserve(nnz);
    ja.reserve(nnz);

    mat.reserve(pat);

    // init matrix, build CSR matrix in parallel
    for (int row=0; row<nrows; ++row) {
        for (int col = -1; col<=1; ++col) {
            int cidx = row + col;
            if (cidx < 0 || cidx >= ncols) continue;

            const double val = (col == 0) ? 2.0 : -1.0;
            values.push_back(val);
            mat.coeffRef(row, cidx) = val;
            ja.push_back(cidx);
        }
    }

    // change matrix
    for (int row=0; row<nrows; ++row) {
        for (int col = -1; col<=1; ++col) {
            int cidx = row + col;
            if (cidx < 0 || cidx >= ncols) continue;

            mat.coeffRef(row, cidx) = 4.0 * mat.coeff(row, cidx);
        }
    }
    // adapt entries of CSR matrix
    for (auto& v : values) v *= 4.0;

    mat.makeCompressed();

    ASSERT_EQ(nnz, mat.nonZeros());
    ASSERT_EQ(nnz, values.size());
    ASSERT_EQ(nnz, ja.size());

    int* ptr = mat.outerIndexPtr();
    int* col = mat.innerIndexPtr();
    double* data = mat.valuePtr();

    for (int r=0; r<(int) nrows; ++r) {
        EXPECT_EQ(ia[r], ptr[r]);
    }

    for (int i=0; i<nnz; ++i) {
        EXPECT_EQ(values[i], data[i]);
        EXPECT_EQ(ja[i], col[i]);
    }
}

