/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file SparseTableCRS.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#ifndef SPARSETABLECRS_H_
#define SPARSETABLECRS_H_

#include <iostream>
#include <vector>
#include <set>
#include <cassert>

#include "Sparsity.h"

namespace MathLib 
{

/**
 * \brief Compressible row storage (CRS) for sparse matrix
 */
template <class T>
struct SparseTableCRS 
{
    /// Dimension of a matrix
    std::size_t dimension;
    /// The number of non-zero entries
    std::size_t nonzero;
    /// Row pointer
    T* row_ptr;
    /// Column index
    T* col_idx;
    /// Data
    double* data;
};

typedef struct SparseTableCRS<unsigned> CRSUnsigned;
typedef struct SparseTableCRS<signed> CRSSigned;

/**
 * convert a row-major sparsity to CRS data
 *
 * \tparam INTTYPE              index type
 * \param row_major_entries     Row-major sparse pattern
 * \return Pointer to newly allocated CRS data
 */
template<class INTTYPE>
MathLib::SparseTableCRS<INTTYPE>* convertRowMajorSparsityToCRS(const RowMajorSparsity &row_major_entries)
{
    const std::size_t n_rows = row_major_entries.size();
    assert(n_rows > 0);
    if (n_rows==0) return 0;
         

    //get number of nonzero 
    std::vector<INTTYPE> *ptr = new std::vector<INTTYPE>(n_rows+1);
    std::vector<INTTYPE> *vec_col_idx = new std::vector<INTTYPE>();
    std::size_t counter_ptr = 0;
    std::size_t cnt_row = 0;

    for (std::size_t i=0; i<n_rows; i++) {
        (*ptr)[cnt_row++] = counter_ptr;         // starting point of the row

        // entries at the i th row
        const std::set<std::size_t> &setConnection = row_major_entries[i];
        //
        for (std::set<std::size_t>::iterator it=setConnection.begin(); it!=setConnection.end(); it++) {
            vec_col_idx->push_back(*it);
            ++counter_ptr;
        }
    }

    (*ptr)[n_rows] = counter_ptr;

    MathLib::SparseTableCRS<INTTYPE> *crs(new MathLib::SparseTableCRS<INTTYPE>);
    crs->dimension = n_rows;
    crs->row_ptr = &(*ptr)[0];
    crs->col_idx = &(*vec_col_idx)[0];
    crs->nonzero = vec_col_idx->size();
    crs->data = new double[vec_col_idx->size()];
    for (std::size_t i=0; i<vec_col_idx->size(); i++)
        crs->data[i] = .0;

    return crs;
}

/**
 * print out CRS data into std::cout
 *
 * \tparam INTTYPE              index type
 * \param crs   CRS data
 */
template<class INTTYPE>
void printoutSparseTableCRS(const MathLib::SparseTableCRS<INTTYPE> *crs)
{
    //output CRS
    std::cout << "PTR:" << std::endl;
    for (std::size_t i=0; i<crs->dimension+1; i++)
        std::cout << crs->row_ptr[i] << ", "; 
    std::cout << std::endl;
    std::cout << "ColID:" << std::endl;
    for (std::size_t i=0; i<crs->nonzero; i++)
        std::cout << crs->col_idx[i] << ", "; 
    std::cout << std::endl;
    std::cout << "Data:" << std::endl;
    for (std::size_t i=0; i<crs->nonzero; i++)
        std::cout << crs->data[i] << ", "; 
    std::cout << std::endl;
}

} // end namespace MathLib

#endif //SPARSETABLECRS_H_

