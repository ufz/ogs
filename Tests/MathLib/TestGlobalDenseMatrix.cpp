/**
 * @file TestGlobalDenseMatrix.cpp
 * @author Thomas Fischer
 * @date Jun 5, 2013
 * @brief 
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <gtest/gtest.h>

#include <limits>

#include "LinAlg/Dense/GlobalDenseMatrix.h"


TEST(MathLib, GlobalDenseMatrix)
{
	const std::size_t n_rows(5);
	const std::size_t n_cols(5);
	MathLib::GlobalDenseMatrix<double, std::size_t> mat0(n_rows,n_cols);
	MathLib::GlobalDenseMatrix<double, std::size_t> mat1(n_rows,n_cols);
	MathLib::GlobalDenseMatrix<double, std::size_t> mat2(n_rows,n_cols-1);

	for (std::size_t i(0); i<n_rows; i++) {
		for (std::size_t j(0); j<n_cols; j++) {
			mat0(i,j) = 1.0 / (i+1.0+j+1.0);
		}
	}

	mat1.setZero();
	mat1 = mat0;
	for (std::size_t i(0); i<n_rows; i++) {
		for (std::size_t j(0); j<n_cols; j++) {
			ASSERT_NEAR(1.0/(i+j+2.0), mat1(i,j), std::numeric_limits<double>::epsilon());
		}
	}

	ASSERT_THROW(mat2 = mat1, std::range_error);

}
