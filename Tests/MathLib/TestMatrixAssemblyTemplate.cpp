/**
 * @file TestMatrixAssemblyTemplate.cpp
 * @author Thomas Fischer
 * @date Jun 13, 2013
 * @brief 
 *
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <gtest/gtest.h>

#include "LinAlg/FinalizeMatrixAssembly.h"
#include "LinAlg/Dense/GlobalDenseMatrix.h"

TEST(MathLib, GlobalDenseMatrixAssembly)
{
	const std::size_t n_rows(5);
	const std::size_t n_cols(5);
	MathLib::GlobalDenseMatrix<double, std::size_t> mat0(n_rows,n_cols);

	// assembly entries
	for (std::size_t i(0); i<n_rows; i++) {
		for (std::size_t j(0); j<n_cols; j++) {
			mat0(i,j) = 1.0 / (i+1.0+j+1.0);
		}
	}

	ASSERT_TRUE(finalizeMatrixAssembly(mat0));
}
