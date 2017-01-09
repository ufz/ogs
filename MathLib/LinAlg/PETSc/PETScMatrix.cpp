/*!
   \file  PETScMatrix.cpp
   \brief Definition of member functions of class PETScMatrix, which provides an interface to
          PETSc matrix routines.

   \author Wenqing Wang
   \date Nov 2011 - Sep 2013

   \copyright
   Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScMatrix.h"

namespace MathLib
{

PETScMatrix::PETScMatrix (const PetscInt nrows, const PETScMatrixOption &mat_opt)
    :_nrows(nrows), _ncols(nrows), _n_loc_rows(PETSC_DECIDE),
     _n_loc_cols(mat_opt.n_local_cols)
{
    if(!mat_opt.is_global_size)
    {
        _n_loc_rows = nrows;
        _n_loc_cols = nrows;
        _nrows = PETSC_DECIDE;
        _ncols = PETSC_DECIDE;
    }

    create(mat_opt.d_nz, mat_opt.o_nz);
}

PETScMatrix::PETScMatrix (const PetscInt nrows, const PetscInt ncols, const PETScMatrixOption &mat_opt)
    :_nrows(nrows), _ncols(ncols),  _n_loc_rows(PETSC_DECIDE),
     _n_loc_cols(mat_opt.n_local_cols)
{
    if(!mat_opt.is_global_size)
    {
        _nrows = PETSC_DECIDE;
        _ncols = PETSC_DECIDE;
        _n_loc_rows = nrows;
        _n_loc_cols = ncols;
    }

    create(mat_opt.d_nz, mat_opt.o_nz);
}

PETScMatrix::PETScMatrix(const PETScMatrix &A)
    : _nrows(A._nrows)
    , _ncols(A._ncols)
    , _n_loc_rows(A._n_loc_rows)
    , _n_loc_cols(A._n_loc_cols)
    , _start_rank(A._start_rank)
    , _end_rank(A._end_rank)
{
    MatConvert(A._A, MATSAME, MAT_INITIAL_MATRIX, &_A);
}

PETScMatrix&
PETScMatrix::operator=(PETScMatrix const& A)
{
    _nrows = A._nrows;
    _ncols = A._ncols;
    _n_loc_rows = A._n_loc_rows;
    _n_loc_cols = A._n_loc_cols;
    _start_rank = A._start_rank;
    _end_rank = A._end_rank;

    if (_A) {
        // TODO this is the slowest option for copying
        MatCopy(A._A, _A, DIFFERENT_NONZERO_PATTERN);
    } else {
        destroy();
        MatConvert(A._A, MATSAME, MAT_INITIAL_MATRIX, &_A);
    }

    return *this;
}

void PETScMatrix::setRowsColumnsZero(std::vector<PetscInt> const& row_pos)
{
    // Each rank (compute core) processes only the rows that belong to the rank itself.
    const PetscScalar one = 1.0;
    const PetscInt nrows = static_cast<PetscInt> (row_pos.size());

    // Each process will only zero its own rows.
    // This avoids all reductions in the zero row routines
    // and thus improves performance for very large process counts.
    // See PETSc doc about MAT_NO_OFF_PROC_ZERO_ROWS.
    MatSetOption(_A, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE);

    // Keep the non-zero pattern for the assignment operator.
    MatSetOption(_A, MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);

    if(nrows>0)
        MatZeroRows(_A, nrows, &row_pos[0], one, PETSC_NULL, PETSC_NULL);
    else
        MatZeroRows(_A, 0, PETSC_NULL, one, PETSC_NULL, PETSC_NULL);
}

void PETScMatrix::viewer(const std::string &file_name, const PetscViewerFormat vw_format)
{
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
    PetscViewerPushFormat(viewer, vw_format);

    finalizeAssembly();

    PetscObjectSetName((PetscObject)_A,"Stiffness_matrix");
    MatView(_A,viewer);

// This preprocessor is only for debugging, e.g. dump the matrix and exit the program.
//#define EXIT_TEST
#ifdef EXIT_TEST
    MatDestroy(_A);
    PetscFinalize();
    exit(0);
#endif

}

void PETScMatrix::create(const PetscInt d_nz, const PetscInt o_nz)
{
    MatCreate(PETSC_COMM_WORLD, &_A);
    MatSetSizes(_A, _n_loc_rows, _n_loc_cols, _nrows, _ncols);

    MatSetFromOptions(_A);

    MatSetType(_A, MATMPIAIJ);
    MatSeqAIJSetPreallocation(_A, d_nz, PETSC_NULL);
    MatMPIAIJSetPreallocation(_A, d_nz, PETSC_NULL, o_nz, PETSC_NULL);
    // If pre-allocation does not work one can use MatSetUp(_A), which is much
    // slower.

    MatGetOwnershipRange(_A, &_start_rank, &_end_rank);
    MatGetSize(_A, &_nrows,  &_ncols);
    MatGetLocalSize(_A, &_n_loc_rows, &_n_loc_cols);
}

bool finalizeMatrixAssembly(PETScMatrix &mat, const MatAssemblyType asm_type)
{
    mat.finalizeAssembly(asm_type);
    return true;
}

} //end of namespace

