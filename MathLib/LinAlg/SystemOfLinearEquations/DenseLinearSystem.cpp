/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-06-25
 * \brief  Definition of the DenseLinearSystem class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DenseLinearSystem.h"

#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"

namespace MathLib
{

/// reset this equation
void DenseLinearSystem::setZero()
{
    const std::size_t dim = getDimension();
    for (std::size_t i=0; i<dim; ++i)
        for (std::size_t j=0; j<dim; ++j)
            _mat(i,j) = .0;
    _rhs = .0;
    _x = .0;
}


/// solve this equation
void DenseLinearSystem::solve()
{
    const std::size_t n_cols = _mat.getNCols();
    for (std::size_t i=0; i<_vec_knownX_id.size(); i++) {
        const std::size_t row_id = _vec_knownX_id[i];
        const double x = _vec_knownX_x[i];
        //A(k, j) = 0.
        for (size_t j=0; j<n_cols; j++)
            _mat(row_id, j) = .0;
        //b_i -= A(i,k)*val, i!=k
        for (size_t j=0; j<n_cols; j++)
            _rhs[j] -= _mat(j, row_id)*x;
        //b_k = val
        _rhs[row_id] = x;
        //A(i, k) = 0., i!=k
        for (size_t j=0; j<n_cols; j++)
            _mat(j, row_id) = .0;
        //A(k, k) = 1.0
        _mat(row_id, row_id) = 1.0; //=x
    }

    MathLib::GaussAlgorithm solver(_mat);
    _x = _rhs;
    solver.execute(&_x[0]);
}

/// printout this equation for debugging
void DenseLinearSystem::printout(std::ostream &os) const
{
    os << "not implemented yet." << std::endl;
}

} // MathLib


