/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Implementation of dense matrix and vector utility functions.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DenseTools.h"

namespace MathLib
{

void applyKnownSolution(DenseMatrix<double> &eqsA, DenseVector<double> &eqsRHS, std::size_t k,
		double val)
{
	const std::size_t n_rows = eqsA.getNRows();
	const std::size_t n_cols = eqsA.getNCols();

	// set all entries of the k-th row of the matrix to zero
	// except the diagonal entry that is set to one
	for (std::size_t j = 0; j < n_cols; j++)
		eqsA(k, j) = .0;
	eqsA(k, k) = 1.0;

	// b_i -= A(i,k)*val, i!=k and
	// set the entries of the k-th column of the matrix to zero
	// except the diagonal entry A(k,k)
	for (std::size_t i = 0; i < k; ++i)
	{
		eqsRHS[i] -= eqsA(i, k) * val;
		eqsA(i, k) = 0.0;
	}
	for (std::size_t i = k + 1; i < n_rows; ++i)
	{
		eqsRHS[i] -= eqsA(i, k) * val;
		eqsA(i, k) = 0.0;
	}

	// b_k = val
	eqsRHS[k] = val;
}

void applyKnownSolution(DenseMatrix<double> &A, DenseVector<double> &b, const std::vector<std::size_t> &vec_knownX_id, const std::vector<double> &vec_knownX_x)
{
	const std::size_t n_bc = vec_knownX_id.size();
	for (std::size_t i_bc = 0; i_bc < n_bc; i_bc++)
	{
		applyKnownSolution(A, b, vec_knownX_id[i_bc], vec_knownX_x[i_bc]);
	}
}


} // MathLib



