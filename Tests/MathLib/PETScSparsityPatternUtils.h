// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <petscmat.h>

#include "MathLib/LinAlg/SparsityPattern.h"

namespace PETScSparsityPatternUtils
{
/// A pattern whose first \c n_filled_rows local rows each hold \c n_block_cols
/// consecutive global columns starting at \c col_start, and whose remaining
/// rows up to \c n_local_rows are empty.
///
/// This is the shape the MPI matrix-interface and linear-solver tests need: a
/// small dense block on the rank's own diagonal sub-matrix, optionally padded
/// with rows that stay structurally empty.
inline MathLib::PETScSparsityPattern denseBlockSparsityPattern(
    PetscInt const n_local_rows, PetscInt const n_filled_rows,
    PetscInt const n_block_cols, PetscInt const col_start)
{
    MathLib::PETScSparsityPattern pattern;
    pattern.row_ptr.reserve(n_local_rows + 1);
    pattern.col_idx.reserve(n_filled_rows * n_block_cols);

    pattern.row_ptr.push_back(0);
    for (PetscInt row = 0; row < n_local_rows; ++row)
    {
        if (row < n_filled_rows)
        {
            for (PetscInt col = 0; col < n_block_cols; ++col)
            {
                pattern.col_idx.push_back(col_start + col);
            }
        }
        pattern.row_ptr.push_back(
            static_cast<PetscInt>(pattern.col_idx.size()));
    }
    return pattern;
}

/// A pattern holding nothing but the diagonal entry of each of the
/// \c n_local_rows local rows, the first of which is global row \c row_start.
inline MathLib::PETScSparsityPattern diagonalSparsityPattern(
    PetscInt const n_local_rows, PetscInt const row_start)
{
    MathLib::PETScSparsityPattern pattern;
    pattern.row_ptr.reserve(n_local_rows + 1);
    pattern.col_idx.reserve(n_local_rows);

    pattern.row_ptr.push_back(0);
    for (PetscInt row = 0; row < n_local_rows; ++row)
    {
        pattern.col_idx.push_back(row_start + row);
        pattern.row_ptr.push_back(
            static_cast<PetscInt>(pattern.col_idx.size()));
    }
    return pattern;
}
}  // namespace PETScSparsityPatternUtils
