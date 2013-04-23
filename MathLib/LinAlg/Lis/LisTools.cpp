/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Implementation of Lis utility functions.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LisTools.h"

#include "logog/include/logog.hpp"

#include "lis.h"

#include "LisMatrix.h"
#include "LisVector.h"

namespace MathLib
{

bool checkLisError(int err)
{
    bool ok = (err == LIS_SUCCESS);
    if (!ok) {
        ERR("***ERROR: Lis error code = %d", err);
    }
    return ok;
}


void applyKnownSolution(LisMatrix &A, LisVector &b, const std::vector<std::size_t> &_vec_knownX_id, const std::vector<double> &_vec_knownX_x)
{
    //Use penalty parameter
    const double penalty_scaling = 1e+10;
    const double _max_diag_coeff = A.getMaxDiagCoeff();
    const double penalty = _max_diag_coeff * penalty_scaling;
    INFO("-> max. absolute value of diagonal entries = %e", _max_diag_coeff);
    INFO("-> penalty scaling = %e", penalty_scaling);
    const std::size_t n_bc = _vec_knownX_id.size();
    for (std::size_t i_bc=0; i_bc<n_bc; i_bc++) {
        const std::size_t rowId = _vec_knownX_id[i_bc];
        const double x = _vec_knownX_x[i_bc];
        //A(k, k) = penalty
        A.setValue(rowId, rowId, penalty);
        //b(k) = x*penalty
        b.set(rowId, x*penalty);
    }
}

} // MathLib



