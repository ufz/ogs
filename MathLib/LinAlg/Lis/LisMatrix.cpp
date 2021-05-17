/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Implementation of the LisMatrix class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LisMatrix.h"

#include <cmath>
#include <cstdlib>

#include "BaseLib/Error.h"
#include "LisCheck.h"
#include "LisVector.h"

namespace MathLib
{
LisMatrix::LisMatrix(std::size_t n_rows, MatrixType mat_type)
    : n_rows_(n_rows),
      mat_type_(mat_type),
      is_assembled_(false),
      use_external_arrays_(false)
{
    int ierr = lis_matrix_create(0, &AA_);
    checkLisError(ierr);
    ierr = lis_matrix_set_size(AA_, 0, n_rows);
    checkLisError(ierr);
    lis_matrix_get_range(AA_, &is_, &ie_);
    ierr = lis_vector_duplicate(AA_, &diag_);
    checkLisError(ierr);
}

LisMatrix::LisMatrix(std::size_t n_rows, int nnz, IndexType* row_ptr,
                     IndexType* col_idx, double* data)
    : n_rows_(n_rows),
      mat_type_(MatrixType::CRS),
      is_assembled_(false),
      use_external_arrays_(true)
{
    int ierr = lis_matrix_create(0, &AA_);
    checkLisError(ierr);
    ierr = lis_matrix_set_size(AA_, 0, n_rows);
    checkLisError(ierr);
    ierr = lis_matrix_set_csr(nnz, row_ptr, col_idx, data, AA_);
    checkLisError(ierr);
    ierr = lis_matrix_assemble(AA_);
    checkLisError(ierr);
    is_assembled_ = true;
    lis_matrix_get_range(AA_, &is_, &ie_);
    ierr = lis_vector_duplicate(AA_, &diag_);
    checkLisError(ierr);
}

LisMatrix::~LisMatrix()
{
    int ierr = 0;
    if (use_external_arrays_)
    {
        ierr = lis_matrix_unset(AA_);
        checkLisError(ierr);
    }
    ierr = lis_matrix_destroy(AA_);
    checkLisError(ierr);
    ierr = lis_vector_destroy(diag_);
    checkLisError(ierr);
}

void LisMatrix::setZero()
{
    // A matrix has to be destroyed and created again because Lis doesn't
    // provide a function to set matrix entries to zero
    int ierr = lis_matrix_destroy(AA_);
    checkLisError(ierr);
    ierr = lis_matrix_create(0, &AA_);
    checkLisError(ierr);
    ierr = lis_matrix_set_size(AA_, 0, n_rows_);
    checkLisError(ierr);
    ierr = lis_vector_set_all(0.0, diag_);
    checkLisError(ierr);

    is_assembled_ = false;
}

int LisMatrix::setValue(IndexType rowId, IndexType colId, double v)
{
    if (v == 0.0)
        return 0;
    lis_matrix_set_value(LIS_INS_VALUE, rowId, colId, v, AA_);
    if (rowId == colId)
        lis_vector_set_value(LIS_INS_VALUE, rowId, v, diag_);
    is_assembled_ = false;
    return 0;
}

int LisMatrix::add(IndexType rowId, IndexType colId, double v)
{
    if (v == 0.0)
        return 0;
    lis_matrix_set_value(LIS_ADD_VALUE, rowId, colId, v, AA_);
    if (rowId == colId)
        lis_vector_set_value(LIS_ADD_VALUE, rowId, v, diag_);
    is_assembled_ = false;
    return 0;
}

void LisMatrix::write(const std::string& filename) const
{
    if (!is_assembled_)
    {
        OGS_FATAL("LisMatrix::write(): matrix not assembled.");
    }
    lis_output_matrix(AA_, LIS_FMT_MM, const_cast<char*>(filename.c_str()));
}

bool finalizeMatrixAssembly(LisMatrix& mat)
{
    LIS_MATRIX& A = mat.getRawMatrix();

    if (!mat.isAssembled())
    {
        int ierr =
            lis_matrix_set_type(A, static_cast<int>(mat.getMatrixType()));
        checkLisError(ierr);
        ierr = lis_matrix_assemble(A);
        checkLisError(ierr);
        mat.is_assembled_ = true;
    }
    return true;
}

}  // namespace MathLib
