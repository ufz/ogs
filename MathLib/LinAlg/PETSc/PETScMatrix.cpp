// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "PETScMatrix.h"

#include "BaseLib/Error.h"
#include "BaseLib/MPI.h"

namespace MathLib
{
PETScMatrix::PETScMatrix(const PetscInt nrows,
                         const PETScSparsityPattern& sparsity_pattern)
    : nrows_(PETSC_DECIDE),
      ncols_(PETSC_DECIDE),
      n_loc_rows_(nrows),
      n_loc_cols_(nrows)
{
    create(sparsity_pattern);
}

PETScMatrix::PETScMatrix(const PetscInt nrows, const PetscInt ncols,
                         const PETScSparsityPattern& sparsity_pattern)
    : nrows_(PETSC_DECIDE),
      ncols_(PETSC_DECIDE),
      n_loc_rows_(nrows),
      n_loc_cols_(ncols)
{
    create(sparsity_pattern);
}

PETScMatrix::PETScMatrix(const PETScMatrix& A)
    : nrows_(A.nrows_),
      ncols_(A.ncols_),
      n_loc_rows_(A.n_loc_rows_),
      n_loc_cols_(A.n_loc_cols_),
      start_rank_(A.start_rank_),
      end_rank_(A.end_rank_)
{
    PetscCallAbort(PETSC_COMM_WORLD,
                   MatConvert(A.A_, MATSAME, MAT_INITIAL_MATRIX, &A_));
}

PETScMatrix& PETScMatrix::operator=(PETScMatrix const& A)
{
    nrows_ = A.nrows_;
    ncols_ = A.ncols_;
    n_loc_rows_ = A.n_loc_rows_;
    n_loc_cols_ = A.n_loc_cols_;
    start_rank_ = A.start_rank_;
    end_rank_ = A.end_rank_;

    if (A_ != nullptr)
    {
        // TODO this is the slowest option for copying
        PetscCallAbort(PETSC_COMM_WORLD,
                       MatCopy(A.A_, A_, DIFFERENT_NONZERO_PATTERN));
    }
    else
    {
        destroy();
        PetscCallAbort(PETSC_COMM_WORLD,
                       MatConvert(A.A_, MATSAME, MAT_INITIAL_MATRIX, &A_));
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
    PetscCallAbort(PETSC_COMM_WORLD,
                   MatSetOption(A_, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE));

    // Keep the non-zero pattern for the assignment operator.
    PetscCallAbort(PETSC_COMM_WORLD,
                   MatSetOption(A_, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE));

    if (nrows > 0)
    {
        PetscCallAbort(PETSC_COMM_WORLD,
                       MatZeroRows(A_, nrows, &row_pos[0], one, PETSC_NULLPTR,
                                   PETSC_NULLPTR));
    }
    else
    {
        PetscCallAbort(PETSC_COMM_WORLD,
                       MatZeroRows(A_, 0, PETSC_NULLPTR, one, PETSC_NULLPTR,
                                   PETSC_NULLPTR));
    }
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
// #define EXIT_TEST
#ifdef EXIT_TEST
    MatDestroy(A_);
    PetscFinalize();
    exit(0);
#endif
}

void PETScMatrix::create(const PETScSparsityPattern& sparsity_pattern)
{
    PetscCallAbort(PETSC_COMM_WORLD, MatCreate(PETSC_COMM_WORLD, &A_));
    PetscCallAbort(PETSC_COMM_WORLD,
                   MatSetSizes(A_, n_loc_rows_, n_loc_cols_, nrows_, ncols_));

    PetscCallAbort(PETSC_COMM_WORLD, MatSetType(A_, MATAIJ));
    PetscCallAbort(PETSC_COMM_WORLD, MatSetFromOptions(A_));

    // Use MATPREALLOCATOR for type-agnostic exact sparsity preallocation.
    // This works regardless of the matrix type chosen via -mat_type.
    preallocateFromSparsityPattern(sparsity_pattern);

    PetscCallAbort(PETSC_COMM_WORLD,
                   MatGetOwnershipRange(A_, &start_rank_, &end_rank_));
    PetscCallAbort(PETSC_COMM_WORLD, MatGetSize(A_, &nrows_, &ncols_));
    PetscCallAbort(PETSC_COMM_WORLD,
                   MatGetLocalSize(A_, &n_loc_rows_, &n_loc_cols_));
}

void PETScMatrix::preallocateFromSparsityPattern(
    const PETScSparsityPattern& sparsity_pattern)
{
    Mat preallocator;
    PetscCallAbort(PETSC_COMM_WORLD,
                   MatCreate(PETSC_COMM_WORLD, &preallocator));
    PetscCallAbort(
        PETSC_COMM_WORLD,
        MatSetSizes(preallocator, n_loc_rows_, n_loc_cols_, nrows_, ncols_));
    PetscCallAbort(PETSC_COMM_WORLD, MatSetType(preallocator, MATPREALLOCATOR));
    PetscCallAbort(PETSC_COMM_WORLD, MatSetUp(preallocator));

    PetscInt row_start;
    PetscInt row_end;
    PetscCallAbort(PETSC_COMM_WORLD,
                   MatGetOwnershipRange(preallocator, &row_start, &row_end));

    PetscInt const n_local = row_end - row_start;

    // The pattern is indexed by local row below, so PETSc's row distribution
    // must agree with the one the pattern was built for. A mismatch would
    // otherwise silently read the wrong rows, or read out of bounds.
    //
    // The condition is rank-local, so it is reduced before aborting: a rank
    // aborting on its own would leave the other ranks hanging in the collective
    // assembly below instead of failing with it.
    if (!BaseLib::MPI::allOf(sparsity_pattern.numberOfRows() == n_local))
    {
        OGS_FATAL(
            "PETScMatrix: the sparsity pattern's local row count and PETSc's "
            "row distribution disagree on at least one rank. On this rank the "
            "pattern describes {} local rows, but PETSc distributed {} rows.",
            sparsity_pattern.numberOfRows(), n_local);
    }

    // Record the nonzero structure of this rank's rows in the MATPREALLOCATOR
    // helper, one row at a time. MATPREALLOCATOR stores positions only, so the
    // values array is ignored and nullptr is passed. All rows here are in
    // [row_start, row_end), so the off-process stash path (which would read
    // the values) is never taken.
    for (PetscInt local_row = 0; local_row < n_local; ++local_row)
    {
        PetscInt const global_row = row_start + local_row;
        PetscInt const ncols = sparsity_pattern.nnzInRow(local_row);
        PetscInt const* const cols = sparsity_pattern.col_idx.data() +
                                     sparsity_pattern.row_ptr[local_row];
        PetscCallAbort(PETSC_COMM_WORLD,
                       MatSetValues(preallocator, 1, &global_row, ncols, cols,
                                    nullptr, INSERT_VALUES));
    }

    PetscCallAbort(PETSC_COMM_WORLD,
                   MatAssemblyBegin(preallocator, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(PETSC_COMM_WORLD,
                   MatAssemblyEnd(preallocator, MAT_FINAL_ASSEMBLY));

    // Insert explicit zeros at every preallocated position by
    // MatPreallocatorPreallocate below (fill = PETSC_TRUE) so that
    // MAT_FINAL_ASSEMBLY in subsequent partial assemblies (e.g. soil-only) does
    // not discard entries that were allocated but not yet written to.
    constexpr PetscBool fill_zeros = PETSC_TRUE;
    PetscCallAbort(PETSC_COMM_WORLD,
                   MatPreallocatorPreallocate(preallocator, fill_zeros, A_));
    // The preallocator is a transient helper, destroyed once the target matrix
    // has its structure. Keeping it would let further matrices built from the
    // same sparsity pattern skip the row-by-row MatSetValues loop above, but
    // the cache would have to live where the pattern lives (the matrix
    // specifications / matrix provider), not in a single matrix. It is not done
    // here: matrices are constructed a handful of times per process, and the
    // loop costs one pass over the pattern, negligible next to assembly.
    PetscCallAbort(PETSC_COMM_WORLD, MatDestroy(&preallocator));
}

bool finalizeMatrixAssembly(PETScMatrix& mat, const MatAssemblyType asm_type)
{
    mat.finalizeAssembly(asm_type);
    return true;
}

void PETScMatrix::addToDiagonal(const PetscScalar value)
{
    // MatShift computes A_ = A_ + value * I, i.e. adds value to every diagonal
    // entry (creating missing ones). Correct in both serial and parallel.
    // Temporarily allow new nonzero entries so a diagonal entry missing from
    // the preallocated sparsity pattern doesn't abort the run.
    PetscCallAbort(
        PETSC_COMM_WORLD,
        MatSetOption(A_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
    PetscCallAbort(PETSC_COMM_WORLD, MatShift(A_, value));
    PetscCallAbort(
        PETSC_COMM_WORLD,
        MatSetOption(A_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
}

void setPreallocationNonzeroOption(PETScMatrix& matrix, bool const has_col_idx)
{
    if (has_col_idx)
    {
        PetscCallAbort(
            PETSC_COMM_WORLD,
            MatSetOption(matrix.getRawMatrix(), MAT_NEW_NONZERO_ALLOCATION_ERR,
                         PETSC_TRUE));
    }
    else
    {
        PetscCallAbort(PETSC_COMM_WORLD,
                       MatSetOption(matrix.getRawMatrix(),
                                    MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE));
    }
}

}  // namespace MathLib
