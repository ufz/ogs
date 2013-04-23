/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Definition of the LisMatrix class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LISMATRIX_H_
#define LISMATRIX_H_

#include <iostream>
#include <cmath>

#include "lis.h"

#include "../MatrixBase.h"

namespace MathLib
{
class LisVector;

/**
 * \brief Lis matrix wrapper class
 */
class LisMatrix : public MatrixBase<double, std::size_t>
{
public:
    /**
     * constructor
     * @param length
     * @param mat_type default 1 (CRS)
     */
    LisMatrix(unsigned length, int mat_type = 1);

    /**
     *
     */
    virtual ~LisMatrix();

    /// reset this equation
    virtual void setZero();

    /// set entry
    virtual int setValue(std::size_t rowId, std::size_t colId, double v)
    {
        if (rowId==colId)
            _max_diag_coeff = std::max(_max_diag_coeff, std::abs(v));
        lis_matrix_set_value(LIS_INS_VALUE, rowId, colId, v, _AA);
        return 0;
    }

    /// add value
    virtual int addValue(std::size_t rowId, std::size_t colId, double v)
    {
        if (rowId==colId)
            _max_diag_coeff = std::max(_max_diag_coeff, std::abs(v));
        lis_matrix_set_value(LIS_ADD_VALUE, rowId, colId, v, _AA);
        return 0;
    }

    /// printout this equation for debugging
    virtual void printout(std::ostream &os=std::cout) const;

    /// get max. diagonal coefficient
    double getMaxDiagCoeff() const { return _max_diag_coeff; };

    /// finish matrix construction
    virtual void finishAssembly();

    /// return if the matrix is already assembled or not
    virtual bool isAssembled() const;

    /// return a raw Lis matrix object
    LIS_MATRIX& getRawMatrix() { return _AA; };

    /// y = mat * x
    void matvec ( const LisVector &x, LisVector &y) const;

private:
    double _max_diag_coeff;
    int _mat_type;
    LIS_MATRIX _AA;
};


} // MathLib

#endif //LISMATRIX_H_

