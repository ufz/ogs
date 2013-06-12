/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Implementation of dense matrix and vector utility functions.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DenseTools.h"

namespace MathLib
{

void applyKnownSolution(DenseMatrix<double> &eqsA, DenseVector<double> &eqsRHS, std::size_t k, double x)
{
    const size_t n_rows = eqsA.getNRows();
    const size_t n_cols = eqsA.getNCols();
    //A(k, j) = 0.
    for (size_t j=0; j<n_cols; j++)
        eqsA(k, j) = .0;
    //A(k, k) = 1,
    eqsA(k,k) = 1.0;
    //b_i -= A(i,k)*val, i!=k
    for (size_t i=0; i<n_rows; ++i)
    {
      if (i==k) continue;
      double v_i_k = .0;
      for (size_t j=0; j<n_cols; ++j)
      {
        if (j==k) {
          v_i_k = eqsA(i, k);
          eqsA(i, k) = .0;
          break;
        } else if (j>k) {
          break;
        }
      }
      eqsRHS[i] -= v_i_k*x;
    }
    //b_k = val
    eqsRHS[k] = x;
}

void applyKnownSolution(DenseMatrix<double> &A, DenseVector<double> &b, const std::vector<std::size_t> &_vec_knownX_id, const std::vector<double> &_vec_knownX_x)
{
    const std::size_t n_bc = _vec_knownX_id.size();
    for (std::size_t i_bc=0; i_bc<n_bc; i_bc++) {
        applyKnownSolution(A, b, _vec_knownX_id[i_bc], _vec_knownX_x[i_bc]);
    }
}


} // MathLib



