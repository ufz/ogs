/*!
   \file
   \brief Definition of member functions of class PETScMatrix, which provides an
   interface to
          PETSc matrix routines.

   \author Wenqing Wang
   \date Nov 2011 - Sep 2013

   \copyright
   Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScMatrix.h"

namespace MathLib
{
PETScMatrix::PETScMatrix(const PetscInt nrows, const PETScMatrixOption& mat_opt)
    : nrows_(nrows),
      ncols_(nrows),
      n_loc_rows_(PETSC_DECIDE),
      n_loc_cols_(mat_opt.n_local_cols)
{
    if (!mat_opt.is_global_size)
    {
        n_loc_rows_ = nrows;
        n_loc_cols_ = nrows;
        nrows_ = PETSC_DECIDE;
        ncols_ = PETSC_DECIDE;
    }

    create(mat_opt.d_nz, mat_opt.o_nz);
}

PETScMatrix::PETScMatrix(const PetscInt nrows, const PetscInt ncols,
                         const PETScMatrixOption& mat_opt)
    : nrows_(nrows),
      ncols_(ncols),
      n_loc_rows_(PETSC_DECIDE),
      n_loc_cols_(mat_opt.n_local_cols)
{
    if (!mat_opt.is_global_size)
    {
        nrows_ = PETSC_DECIDE;
        ncols_ = PETSC_DECIDE;
        n_loc_rows_ = nrows;
        n_loc_cols_ = ncols;
    }

    create(mat_opt.d_nz, mat_opt.o_nz);
}

PETScMatrix::PETScMatrix(const PETScMatrix& A)
    : nrows_(A.nrows_),
      ncols_(A.ncols_),
      n_loc_rows_(A.n_loc_rows_),
      n_loc_cols_(A.n_loc_cols_),
      start_rank_(A.start_rank_),
      end_rank_(A.end_rank_)
{
    MatConvert(A.A_, MATSAME, MAT_INITIAL_MATRIX, &A_);
}

PETScMatrix& PETScMatrix::operator=(PETScMatrix const& A)
{
    nrows_ = A.nrows_;
    ncols_ = A.ncols_;
    n_loc_rows_ = A.n_loc_rows_;
    n_loc_cols_ = A.n_loc_cols_;
    start_rank_ = A.start_rank_;
    end_rank_ = A.end_rank_;

    if (A_)
    {
        // TODO this is the slowest option for copying
        MatCopy(A.A_, A_, DIFFERENT_NONZERO_PATTERN);
    }
    else
    {
        destroy();
        MatConvert(A.A_, MATSAME, MAT_INITIAL_MATRIX, &A_);
    }

    return *this;
}

void PETScMatrix::setRowsColumnsZero(std::vector<PetscInt> const& row_pos)
{
    // Each rank (compute core) processes only the rows that belong to the rank
    // itself.
    const PetscScalar one = 1.0;
    const PetscInt nrows = static_cast<PetscInt>(row_pos.size());

    // Each process will only zero its own rows.
    // This avoids all reductions in the zero row routines
    // and thus improves performance for very large process counts.
    // See PETSc doc about MAT_NO_OFF_PROC_ZERO_ROWS.
    MatSetOption(A_, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE);

    // Keep the non-zero pattern for the assignment operator.
    MatSetOption(A_, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);

    if (nrows > 0)
        MatZeroRows(A_, nrows, &row_pos[0], one, PETSC_NULL, PETSC_NULL);
    else
        MatZeroRows(A_, 0, PETSC_NULL, one, PETSC_NULL, PETSC_NULL);
}

void PETScMatrix::viewer(const std::string& file_name,
                         const PetscViewerFormat vw_format)
{
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
    PetscViewerPushFormat(viewer, vw_format);

    finalizeAssembly();

    PetscObjectSetName((PetscObject)A_, "Stiffness_matrix");
    MatView(A_, viewer);

// This preprocessor is only for debugging, e.g. dump the matrix and exit the
// program.
//#define EXIT_TEST
#ifdef EXIT_TEST
    MatDestroy(A_);
    PetscFinalize();
    exit(0);
#endif
}

void PETScMatrix::create(const PetscInt d_nz, const PetscInt o_nz)
{
    MatCreate(PETSC_COMM_WORLD, &A_);
    MatSetSizes(A_, n_loc_rows_, n_loc_cols_, nrows_, ncols_);

    MatSetType(A_, MATAIJ);
    MatSetFromOptions(A_);

    MatSeqAIJSetPreallocation(A_, d_nz, PETSC_NULL);
    MatMPIAIJSetPreallocation(A_, d_nz, PETSC_NULL, o_nz, PETSC_NULL);
    // If pre-allocation does not work one can use MatSetUp(A_), which is much
    // slower.

    MatGetOwnershipRange(A_, &start_rank_, &end_rank_);
    MatGetSize(A_, &nrows_, &ncols_);
    MatGetLocalSize(A_, &n_loc_rows_, &n_loc_cols_);
}

bool finalizeMatrixAssembly(PETScMatrix& mat, const MatAssemblyType asm_type)
{
    mat.finalizeAssembly(asm_type);
    return true;
}

}  // namespace MathLib
