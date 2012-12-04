/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file Sparsity.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#ifndef SPARSITY_H_
#define SPARSITY_H_

#include <vector>
#include <set>

namespace MathLib 
{

/**
 * \brief Row-major sparse pattern
 */
typedef std::vector<std::set<std::size_t> > RowMajorSparsity;

/**
 * convert a row-major sparsity to CRS data
 *
 * @tparam INTTYPE              index type
 * @param row_major_entries     Row-major sparse pattern
 * @param n_rows                The number of rows
 * @param row_ptr               Pointer to row index
 * @param col_idx               Pointer to column index
 * @param nonzero               The number of non-zero entries
 * @param data                  Pointer to non-zero values
 */
template<class INTTYPE>
void convertRowMajorSparsityToCRS(const RowMajorSparsity &row_major_entries, std::size_t &n_rows, INTTYPE* &row_ptr, INTTYPE* &col_idx, std::size_t &nonzero, double* &data);

} // end namespace MathLib

#include "Sparsity.tpp"

#endif // SPARSITY_H_
