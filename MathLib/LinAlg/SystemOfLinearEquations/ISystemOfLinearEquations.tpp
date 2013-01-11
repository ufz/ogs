/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-06-25
 * \brief  Implementation of the ISystemOfLinearEquations class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ISYSTEMOFLINEAREQUATIONS_TPP_
#define ISYSTEMOFLINEAREQUATIONS_TPP_

namespace MathLib
{

template <class T_DENSE_MATRIX>
void ISystemOfLinearEquations::addSubMat(const std::vector<std::size_t> &vec_row_pos, const std::vector<std::size_t> &vec_col_pos, const T_DENSE_MATRIX &sub_matrix, double fkt)
{
    const std::size_t n_rows = vec_row_pos.size();
    const std::size_t n_cols = vec_col_pos.size();
    for (std::size_t i=0; i<n_rows; i++) {
        const std::size_t rowId = vec_row_pos[i];
        if (rowId==index_npos) continue;
        for (std::size_t j=0; j<n_cols; j++) {
            const std::size_t colId = vec_col_pos[j];
            if (colId==index_npos) continue;
            addMatEntry(rowId, colId, fkt*sub_matrix(i,j));
        }
    }
}

template <class T_DENSE_MATRIX>
void ISystemOfLinearEquations::addSubMat(const std::vector<std::size_t> &vec_pos, const T_DENSE_MATRIX &sub_matrix, double fkt)
{
    addSubMat(vec_pos, vec_pos, sub_matrix, fkt);
}

inline void ISystemOfLinearEquations::addSubRHS(const std::vector<std::size_t> &vec_row_pos, const double *sub_vector, double fkt)
{
    for (std::size_t i=0; i<vec_row_pos.size(); i++) {
        const std::size_t rowId = vec_row_pos[i];
        if (rowId==index_npos) continue;
        addRHSVec(rowId, sub_vector[i]*fkt);
    }
}

template <class T_DENSE_VECTOR>
void ISystemOfLinearEquations::addSubRHS(const std::vector<std::size_t> &vec_row_pos, const T_DENSE_VECTOR &sub_vector, double fkt)
{
    for (std::size_t i=0; i<vec_row_pos.size(); i++) {
        const std::size_t rowId = vec_row_pos[i];
        if (rowId==index_npos) continue;
        addRHSVec(rowId, sub_vector(i)*fkt);
    }
}

}

#endif //ISYSTEMOFLINEAREQUATIONS_TPP_

