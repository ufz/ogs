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

LisMatrix::LisMatrix(std::size_t dimension, LisOption::MatrixType mat_type)
: _n_rows(dimension), _max_diag_coeff(.0), _mat_type(mat_type)
{
    int ierr = 0;
    ierr = lis_matrix_create(0, &_AA); checkLisError(ierr);
    ierr = lis_matrix_set_size(_AA, 0, dimension); checkLisError(ierr);
    lis_matrix_get_range(_AA, &_is, &_ie);
}

LisMatrix::~LisMatrix()
{
    int ierr = 0;
    ierr = lis_matrix_destroy(_AA); checkLisError(ierr);
}

void LisMatrix::setZero()
{
    int ierr = 0;
    // A matrix has to be destroyed and created again because Lis doesn't provide a
    // function to set matrix entries to zero
    ierr = lis_matrix_destroy(_AA); checkLisError(ierr);
    ierr = lis_matrix_create(0, &_AA); checkLisError(ierr);
    ierr = lis_matrix_set_size(_AA, 0, _n_rows); checkLisError(ierr);

    _max_diag_coeff = .0;
}

void LisMatrix::write(const std::string &filename) const
{
    lis_output_matrix(_AA, LIS_FMT_MM, const_cast<char*>(filename.c_str()));
}

void LisMatrix::matvec (const LisVector &x, LisVector &y) const
{
    int ierr = lis_matvec(_AA, const_cast<LisVector*>(&x)->getRawVector(), y.getRawVector()); checkLisError(ierr);
}

void finishMatrixAssembly(LisMatrix &mat)
{
    int ierr = lis_matrix_set_type(mat._AA, static_cast<int>(mat._mat_type)); checkLisError(ierr);
    ierr = lis_matrix_assemble(mat._AA); checkLisError(ierr);
};

bool isMatrixAssembled(LisMatrix &mat)
{
	return lis_matrix_is_assembled(mat._AA);
}

} //MathLib
