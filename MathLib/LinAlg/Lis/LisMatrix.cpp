/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Implementation of the LisMatrix class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LisMatrix.h"

#include "LisVector.h"
#include "LisCheck.h"

namespace MathLib
{

LisMatrix::LisMatrix(std::size_t n_rows, LisOption::MatrixType mat_type)
	: _n_rows(n_rows), _max_diag_coeff(.0), _mat_type(mat_type), _is_assembled(false)
{
    int ierr = lis_matrix_create(0, &_AA);
    checkLisError(ierr);
    ierr = lis_matrix_set_size(_AA, 0, n_rows);
    checkLisError(ierr);
    lis_matrix_get_range(_AA, &_is, &_ie);
}

LisMatrix::~LisMatrix()
{
    int ierr = lis_matrix_destroy(_AA);
    checkLisError(ierr);
}

void LisMatrix::setZero()
{
    // A matrix has to be destroyed and created again because Lis doesn't provide a
    // function to set matrix entries to zero
    int ierr = lis_matrix_destroy(_AA);
    checkLisError(ierr);
    ierr = lis_matrix_create(0, &_AA);
    checkLisError(ierr);
    ierr = lis_matrix_set_size(_AA, 0, _n_rows);
    checkLisError(ierr);

    _max_diag_coeff = .0;
}

int LisMatrix::setValue(std::size_t rowId, std::size_t colId, double v)
{
    if (rowId==colId)
        _max_diag_coeff = std::max(_max_diag_coeff, std::abs(v));
    lis_matrix_set_value(LIS_INS_VALUE, rowId, colId, v, _AA);
    return 0;
}

int LisMatrix::addValue(std::size_t rowId, std::size_t colId, double v)
{
    if (rowId==colId)
        _max_diag_coeff = std::max(_max_diag_coeff, std::abs(v));
    lis_matrix_set_value(LIS_ADD_VALUE, rowId, colId, v, _AA);
    return 0;
}

void LisMatrix::write(const std::string &filename) const
{
    lis_output_matrix(_AA, LIS_FMT_MM, const_cast<char*>(filename.c_str()));
}

void LisMatrix::matvec (const LisVector &x, LisVector &y) const
{
    int ierr = lis_matvec(_AA, const_cast<LisVector*>(&x)->getRawVector(), y.getRawVector());
    checkLisError(ierr);
}

bool finalizeMatrixAssembly(LisMatrix &mat)
{
    LIS_MATRIX &A = mat.getRawMatrix();
    // commented out below because lis_matrix_is_assembled() always returns the same value.
    //    if (LIS_SUCCESS!=lis_matrix_is_assembled(A)) {
    if (!mat.isAssembled()) {
        int ierr = lis_matrix_set_type(A, static_cast<int>(mat.getMatrixType()));
        checkLisError(ierr);
        ierr = lis_matrix_assemble(A); //checkLisError(ierr);
        mat._is_assembled = true;
   }
    return true;
};


} //MathLib
