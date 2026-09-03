// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

#ifdef USE_PETSC
#include <petscsystypes.h>
#endif

namespace MathLib
{
/// Per-row nonzero count pattern for the Eigen (serial) build.
///
/// \c number_non_zeros_per_row is sized to the number of global matrix rows
/// and holds the number of nonzeros per row.
template <typename IndexType>
struct SerialSparsityPattern
{
    std::vector<IndexType> number_non_zeros_per_row;

    /// Number of nonzeros in the given global row.
    IndexType nnzInRow(IndexType const row) const
    {
        return number_non_zeros_per_row[row];
    }

    /// Number of global rows the pattern describes.
    IndexType numberOfRows() const
    {
        return static_cast<IndexType>(number_non_zeros_per_row.size());
    }
};

#ifdef USE_PETSC
/// Sparsity pattern for the PETSc (parallel) build.
///
/// CSR representation with global column indices sorted within each local row.
/// Used for type-agnostic preallocation via MATPREALLOCATOR. The index type is
/// PetscInt, because that is the only type PETSc accepts.
struct PETScSparsityPattern
{
    /// CSR row pointers (length n_local_rows + 1).
    std::vector<PetscInt> row_ptr;
    /// Global column indices, sorted within each local row.
    std::vector<PetscInt> col_idx;

    /// Number of nonzeros in the given local row.
    PetscInt nnzInRow(PetscInt const row) const
    {
        return row_ptr[row + 1] - row_ptr[row];
    }

    /// Number of local rows the pattern describes.
    PetscInt numberOfRows() const
    {
        // A default-constructed pattern has an empty row_ptr and describes
        // zero rows; only a populated pattern carries the trailing sentinel.
        return row_ptr.empty() ? 0 : static_cast<PetscInt>(row_ptr.size()) - 1;
    }
};
#endif
}  // namespace MathLib
