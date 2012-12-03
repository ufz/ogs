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

    //get number of nonzero
    row_ptr = new INTTYPE[n_rows+1];
    std::vector<INTTYPE> vec_col_idx;
    size_t counter_ptr = 0;
    size_t cnt_row = 0;

    for (size_t i=0; i<n_rows; i++) {
        row_ptr[cnt_row++] = counter_ptr;         // starting point of the row

        // entries at the i th row
        const std::set<size_t> &setConnection = row_major_entries[i];
        //
        for (std::set<size_t>::iterator it=setConnection.begin(); it!=setConnection.end(); it++) {
            vec_col_idx.push_back(*it);
            ++counter_ptr;
        }
    }

    row_ptr[n_rows] = counter_ptr;
    nonzero = vec_col_idx.size();
    col_idx = new INTTYPE[vec_col_idx.size()];
    for (size_t i=0; i<nonzero; i++)
        col_idx[i] = vec_col_idx[i];
    data = new double[nonzero];
    for (size_t i=0; i<nonzero; i++)
        data[i] = .0;
}

} // end namespace MathLib

#endif // SPARSITY_TPP_
