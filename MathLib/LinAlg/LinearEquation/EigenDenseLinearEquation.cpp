/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file EigenDenseLinearEquation.cpp
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#include "EigenDenseLinearEquation.h"


namespace MathLib
{

EigenDenseLinearEquation::EigenDenseLinearEquation(std::size_t length, RowMajorSparsity* sp)
: ISystemOfLinearEquations(length, sp)
{
    resize(length);
}

void EigenDenseLinearEquation::resize(std::size_t length)
{
    _A = MatrixType::Zero(length, length);
    _b = VectorType::Zero(length);
    _x = VectorType::Zero(length);
    //reset();
}

void EigenDenseLinearEquation::setZero()
{
    _A *= .0;
    _b *= .0;
    _x *= .0;
}

void EigenDenseLinearEquation::setKnownSolution(std::size_t row_id, double x)
{
    const std::size_t n_cols = _A.cols();
    //A(k, j) = 0.
    for (std::size_t j=0; j<n_cols; j++)
        _A(row_id, j) = .0;
    //b_i -= A(i,k)*val, i!=k
    for (std::size_t j=0; j<n_cols; j++)
        _b[j] -= _A(j, row_id)*x;
    //b_k = val
    _b[row_id] = x;
    //A(i, k) = 0., i!=k
    for (std::size_t j=0; j<n_cols; j++)
        _A(j, row_id) = .0;
    //A(k, k) = 1.0
    _A(row_id, row_id) = 1.0; //=x
}

void EigenDenseLinearEquation::setKnownSolution(const std::vector<std::size_t> &vec_id, const std::vector<double> &vec_x)
{
    for (std::size_t i=0; i<vec_id.size(); ++i)
        setKnownSolution(vec_id[i], vec_x[i]);
}

void EigenDenseLinearEquation::printout(std::ostream &os) const
{
    os << "A=" << _A << std::endl;
    os << "x=" << _x << std::endl;
    os << "b=" << _b << std::endl;
}

void EigenDenseLinearEquation::solve()
{
    _x = _A.colPivHouseholderQr().solve(_b);
}

}
