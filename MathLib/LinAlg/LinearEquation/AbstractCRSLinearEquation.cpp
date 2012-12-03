/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file AbstractCRSLinearEquation.cpp
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#include "AbstractCRSLinearEquation.h"

#include <algorithm>

namespace MathLib
{
template class AbstractCRSLinearEquation<signed>;
template class AbstractCRSLinearEquation<unsigned>;

template<typename IDX_TYPE> void AbstractCRSLinearEquation<IDX_TYPE>::create(size_t length, RowMajorSparsity *sparsity)
{
    std::size_t dim = 0;
    std::size_t nonzero = 0;
    IDX_TYPE* row_ptr = nullptr;
    IDX_TYPE* col_idx = nullptr;
    double* data = nullptr;
    convertRowMajorSparsityToCRS<IDX_TYPE>(*sparsity, dim, row_ptr, col_idx, nonzero, data);
    assert (length == dim);
    _A = new CRSMatrix<double, IDX_TYPE>(dim, row_ptr, col_idx, data);
    _b.resize(length);
    _x.resize(length);
}

template<typename IDX_TYPE> void AbstractCRSLinearEquation<IDX_TYPE>::setZero()
{
    (*_A) = .0;
    _b.assign(_b.size(), .0);
    _x.assign(_x.size(), .0);
    _vec_knownX_id.clear();
    _vec_knownX_x.clear();
}


template<typename IDX_TYPE> void AbstractCRSLinearEquation<IDX_TYPE>::solve()
{
    if (_vec_knownX_id.size()>0) {
        CRSMatrix<double, IDX_TYPE>* tmp_A = new CRSMatrix<double, IDX_TYPE>(*getMatEntry());
        double *org_eqsRHS = getRHSVec();
        double *org_eqsX = getSolVec();
        std::vector<double> _tmp_b;
        std::vector<double> _tmp_x;
        std::map<size_t,size_t> _map_solved_orgEqs;

        //std::cout << "#before\n";
        //this->printout();
        setKnownXi_ReduceSizeOfEQS(tmp_A, org_eqsRHS, org_eqsX, _vec_knownX_id, _vec_knownX_x, _tmp_b, _tmp_x, _map_solved_orgEqs);
        //std::cout << "\n#after\n";
        //tmp_A->printMat();

        solveEqs(tmp_A, &_tmp_b[0], &_tmp_x[0]);

        const size_t dim = tmp_A->getNRows();
        for (size_t i=0; i<dim; i++) {
            setSolVec(_map_solved_orgEqs[i], _tmp_x[i]);
        }

        delete tmp_A;
    } else {
        solveEqs(getMatEntry(), getRHSVec(), getSolVec());
    }
}

template<typename IDX_TYPE> void AbstractCRSLinearEquation<IDX_TYPE>::setKnownXi_ReduceSizeOfEQS(CRSMatrix<double, IDX_TYPE> *A, double *org_eqsRHS, double *org_eqsX, const std::vector<size_t> &vec_id, const std::vector<double> &vec_x, std::vector<double> &out_b, std::vector<double> &out_x, std::map<size_t,size_t> &map_solved_orgEqs)
{
    assert(vec_id.size()==vec_x.size());

    const size_t n_org_rows = A->getNRows();

    std::vector<IDX_TYPE> removed_rows(vec_id.size());
    for (size_t i=0; i<vec_id.size(); i++) {
        const size_t id = vec_id[i];
        const double val = vec_x[i];
        removed_rows[i] = id;

        //b_i -= A(i,k)*val, i!=k
        for (IDX_TYPE j=0; j<A->getNCols(); j++)
            org_eqsRHS[j] -= A->getValue(j, id)*val;
        //b_k = A_kk*val
        org_eqsRHS[id] = val; //=eqsA(id, id)*val;
        org_eqsX[id] = val; //=eqsA(id, id)*val;
    }

    //remove rows and columns
    std::sort(removed_rows.begin(), removed_rows.end());
    A->eraseEntries(removed_rows.size(), &removed_rows[0]);
    const size_t n_new_rows = n_org_rows-removed_rows.size();

    //remove X,RHS
    out_b.resize(n_new_rows);
    out_x.resize(n_new_rows);
    size_t new_id = 0;
    for (size_t i=0; i<n_org_rows; i++) {
        if (std::find(removed_rows.begin(), removed_rows.end(), i)!=removed_rows.end()) continue;
        out_b[new_id] = org_eqsRHS[i];
        out_x[new_id] = org_eqsX[i];
        map_solved_orgEqs[new_id] = i;
        new_id++;
    }
}
} //end
