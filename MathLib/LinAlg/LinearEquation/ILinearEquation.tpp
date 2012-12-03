/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file ILinearEquation.cpp
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#include "ILinearEquation.h"

namespace MathLib
{

void ILinearEquation::addAsub(const std::vector<size_t> &vec_row_pos, const std::vector<size_t> &vec_col_pos, const LocalMatrix &sub_matrix, double fkt)
{
    const size_t n_rows = vec_row_pos.size();
    const size_t n_cols = vec_col_pos.size();
    for (size_t i=0; i<n_rows; i++) {
        const size_t rowId = vec_row_pos[i];
        if (rowId==index_npos) continue;
        for (size_t j=0; j<n_cols; j++) {
            const size_t colId = vec_col_pos[j];
            if (colId==index_npos) continue;
            addA(rowId, colId, fkt*sub_matrix(i,j));
        }
    }
}

void ILinearEquation::addAsub(const std::vector<size_t> &vec_pos, const LocalMatrix &sub_matrix, double fkt)
{
    addAsub(vec_pos, vec_pos, sub_matrix, fkt);
}

void ILinearEquation::addRHSsub(const std::vector<size_t> &vec_row_pos, const double *sub_vector, double fkt)
{
    for (size_t i=0; i<vec_row_pos.size(); i++) {
        const size_t rowId = vec_row_pos[i];
        if (rowId==index_npos) continue;
        addRHS(rowId, sub_vector[i]*fkt);
    }
}

void ILinearEquation::addRHSsub(const std::vector<size_t> &vec_row_pos, const LocalVector &sub_vector, double fkt)
{
    for (size_t i=0; i<vec_row_pos.size(); i++) {
        const size_t rowId = vec_row_pos[i];
        if (rowId==index_npos) continue;
        addRHS(rowId, sub_vector(i)*fkt);
    }
}

}
