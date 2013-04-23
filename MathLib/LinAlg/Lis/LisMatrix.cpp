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
#include "LisTools.h"

namespace MathLib
{

LisMatrix::LisMatrix(unsigned dimension, int mat_type)
: MatrixBase<double, std::size_t>(dimension), _max_diag_coeff(.0), _mat_type(mat_type)
{
    int ierr = 0;
    ierr = lis_matrix_create(0, &_AA); checkLisError(ierr);
    ierr = lis_matrix_set_size(_AA, 0, dimension); checkLisError(ierr);
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

    _max_diag_coeff = .0;
}

bool LisMatrix::isAssembled() const
{
    return lis_matrix_is_assembled(_AA);
}

void LisMatrix::finishAssembly()
{
    int ierr = lis_matrix_set_type(_AA, _mat_type); checkLisError(ierr);
    ierr = lis_matrix_assemble(_AA); checkLisError(ierr);
}

void LisMatrix::printout(std::ostream &os) const
{
    os << "#A=" << std::endl;
}

void LisMatrix::matvec (const LisVector &x, LisVector &y) const
{
    int ierr = lis_matvec(_AA, const_cast<LisVector*>(&x)->getRawVector(), y.getRawVector()); checkLisError(ierr);
}

} //MathLib
