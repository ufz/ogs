/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file Sparsity.tpp
 *
 * Created on 2012-12-03 by Norihiro Watanabe
 */

#ifndef SPARSITY_TPP_
#define SPARSITY_TPP_

#include <cassert>

namespace MathLib 
{

template<class INTTYPE>
void convertRowMajorSparsityToCRS(const RowMajorSparsity &row_major_entries, std::size_t &n_rows, INTTYPE* &row_ptr, INTTYPE* &col_idx, std::size_t &nonzero, double* &data)
{
    n_rows = row_major_entries.size();
    assert(n_rows > 0);
    if (n_rows==0) return;

    // get the number of nonzero entries and store the locations in the nonzero
    // vector that start a row
    nonzero = 0;
    row_ptr = new INTTYPE[n_rows+1];
    for (std::size_t i=0; i<n_rows; i++) {
        row_ptr[i] = nonzero;         // starting point of the row
        nonzero += row_major_entries[i].size(); // entries at the i th row
    }
    row_ptr[n_rows] = nonzero;

    // store column indexes of nonzero entries
    col_idx = new INTTYPE[nonzero];
    size_t cnt_entries = 0;
    for (std::size_t i=0; i<n_rows; i++) {
        const std::set<std::size_t> &setConnection = row_major_entries[i];
        for (std::set<std::size_t>::iterator it=setConnection.begin(); it!=setConnection.end(); ++it) {
            col_idx[cnt_entries++] = *it;
        }
    }

    // allocate memory for nonzero entries
    data = new double[nonzero];
    for (std::size_t i=0; i<nonzero; i++)
        data[i] = .0;
}

} // end namespace MathLib

#endif // SPARSITY_TPP_
