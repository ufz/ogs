/*!
   \file
   \brief Declaration of class PETScMatrix, which provides an interface to
          PETSc matrix routines.

   \author Wenqing Wang
   \date Nov 2013 - 2014

   \copyright
    Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#pragma once

#include <string>
#include <vector>

#include "MathLib/LinAlg/RowColumnIndices.h"
#include "PETScMatrixOption.h"
#include "PETScVector.h"

typedef Mat PETSc_Mat;

namespace MathLib
{
/*!
   \brief Wrapper class for PETSc matrix routines for matrix.
*/
class PETScMatrix
{
public:
    using IndexType = PetscInt;

public:
    PETScMatrix() {}
    /*!
      \brief        Constructor for a square matrix partitioning with more
      options
      \param nrows  The number of rows of the matrix or the local matrix.
      \param mat_op The configuration information for creating a matrix.
    */
    PETScMatrix(const PetscInt nrows,
                const PETScMatrixOption& mat_op = PETScMatrixOption());

    /*!
      \brief        Constructor for a rectangular matrix partitioning with more
      options
      \param nrows  The number of global or local rows.
      \param ncols  The number of global or local columns.
      \param mat_op The configuration information for creating a matrix.
    */
    PETScMatrix(const PetscInt nrows, const PetscInt ncols,
                const PETScMatrixOption& mat_op = PETScMatrixOption());

    ~PETScMatrix() { destroy(); }
    PETScMatrix(PETScMatrix const& A);

    PETScMatrix& operator=(PETScMatrix const& A);

    /*!
       \brief          Perform MPI collection of assembled entries in buffer
       \param asm_type Assembly type, either MAT_FLUSH_ASSEMBLY
                       or MAT_FINAL_ASSEMBLY
    */
    void finalizeAssembly(const MatAssemblyType asm_type = MAT_FINAL_ASSEMBLY)
    {
        MatAssemblyBegin(A_, asm_type);
        MatAssemblyEnd(A_, asm_type);
    }

    /// Get the number of rows.
    PetscInt getNumberOfRows() const { return nrows_; }
    /// Get the number of columns.
    PetscInt getNumberOfColumns() const { return ncols_; }
    /// Get the number of local rows.
    PetscInt getNumberOfLocalRows() const { return n_loc_rows_; }
    /// Get the number of local columns.
    PetscInt getNumberOfLocalColumns() const { return n_loc_cols_; }
    /// Get the start global index of the rows of the same rank.
    PetscInt getRangeBegin() const { return start_rank_; }
    /// Get the end global index of the rows in the same rank.
    PetscInt getRangeEnd() const { return end_rank_; }
    /// Get matrix reference.
    Mat& getRawMatrix() { return A_; }
    /*! Get a matrix reference.
     *
     * \warning
     * This method is dangerous insofar as you can do arbitrary things also
     * with a const PETSc matrix.
     */
    Mat const& getRawMatrix() const { return A_; }
    /// Set all entries to zero.
    void setZero() { MatZeroEntries(A_); }
    /*!
       \brief Set the specified rows to zero except diagonal entries, i.e.
              \f$A(k, j) = \begin{cases}
                0.0, &j\not=k, j=1,2,\dots,k-1, k+1, \dots, n \\
                1.0, &j = k
              \end{cases}\f$, where \f$k \in \mbox{row\_pos}\f$
              This function must be called by all ranks.
       \param row_pos The row indices of the specified rows.
    */
    void setRowsColumnsZero(std::vector<PetscInt> const& row_pos);

    /*!
       \brief       Set a single entry with a value.
       \param i     The row index.
       \param j     The column index.
       \param value The entry value.
    */
    void set(const PetscInt i, const PetscInt j, const PetscScalar value)
    {
        MatSetValue(A_, i, j, value, INSERT_VALUES);
    }

    /*!
       \brief Add value to a single entry.
       \param i     The row index.
       \param j     The column index.
       \param value The entry value.
    */
    void add(const PetscInt i, const PetscInt j, const PetscScalar value)
    {
        MatSetValue(A_, i, j, value, ADD_VALUES);
    }

    /*!
       \brief Add sub-matrix at positions given by global \c indices, in which
       negative index indicates ghost entry.

       In order to use MatZeroRows to apply Dirichlet boundary condition,
       entries in the rows with the negative global indices are skipped to added
       to the global matrix, meanwhile entries in the columns with the negative
       global indices are added the global matrix. By using MatZeroRows to apply
       Dirichlet boundary condition, the off diagonal entries in ghost rows of
       the global matrix are set to zero, while the off diagonal entries in
       ghost rows of the global matrix are assembled and kept for linear solver.

       For the setting of Dirichlet boundary condition
       in PETSc, please refer to the
       [PETSc:Documentation:FAQ](https://petsc.org/release/faq/#when-solving-a-system-with-dirichlet-boundary-conditions-i-can-use-matzerorows-to-eliminate-the-dirichlet-rows-but-this-results-in-a-non-symmetric-system-how-can-i-apply-dirichlet-boundary-conditions-but-keep-the-matrix-symmetric).
     */
    template <class T_DENSE_MATRIX>
    void add(RowColumnIndices<PetscInt> const& indices,
             const T_DENSE_MATRIX& sub_matrix)
    {
        // Set global column indices to positive to allow all entries of columns
        // to be added to the global matrix. For the ghost columns, only the
        // off diagonal entries are added due to the negative indices of the
        // corresponding rows.
        std::vector<PetscInt> cols;
        cols.reserve(indices.columns.size());
        for (auto col : indices.columns)
        {
            // Ghost entries, and its original index is 0.
            if (col == -ncols_)
            {
                cols.push_back(0);
            }
            else
            {
                cols.push_back(std::abs(col));
            }
        }

        add(indices.rows, cols, sub_matrix);
    }

    /*!
        \brief          Add a dense sub-matrix to a PETSc matrix.
        \param row_pos  The global indices of the rows of the dense sub-matrix.
        \param col_pos  The global indices of the columns of the dense
                        sub-matrix.
       \param sub_mat  A dense sub-matrix to be added.
    */
    template <class T_DENSE_MATRIX>
    void add(std::vector<PetscInt> const& row_pos,
             std::vector<PetscInt> const& col_pos,
             const T_DENSE_MATRIX& sub_mat);

    /*! View the global vector for test purpose. Do not use it for output a big
       vector.
        \param file_name  File name for output
        \param vw_format  File format listed as:
         PETSC_VIEWER_DEFAULT            Default format
         PETSC_VIEWER_ASCII_MATLAB       MATLAB format
         PETSC_VIEWER_ASCII_DENSE        Print matrix as dense
         PETSC_VIEWER_ASCII_IMPL         Implementation-specific format
                                         (which is in many cases the same as
                                         the default)
         PETSC_VIEWER_ASCII_INFO         Basic information about object
         PETSC_VIEWER_ASCII_INFO_DETAIL  More detailed info about object
         PETSC_VIEWER_ASCII_COMMON       Identical output format for all objects
                                         of a particular type
         PETSC_VIEWER_ASCII_INDEX        (for vectors) Prints the vector element
                                         number next to each vector entry
         PETSC_VIEWER_ASCII_SYMMODU      Print parallel vectors without
                                         indicating the processor ranges
         PETSC_VIEWER_ASCII_VTK          Outputs the object to a VTK file
         PETSC_VIEWER_NATIVE             Store the object to the binary file in
                                         its native format (for example,
                                         dense matrices are stored as dense),
                                         DMDA vectors are dumped directly to
                                         the file instead of being first put in
                                         the natural ordering
         PETSC_VIEWER_DRAW_BASIC         Views the vector with a simple 1d plot
         PETSC_VIEWER_DRAW_LG            Views the vector with a line graph
         PETSC_VIEWER_DRAW_CONTOUR       Views the vector with a contour plot
    */
    void viewer(const std::string& file_name,
                const PetscViewerFormat vw_format = PETSC_VIEWER_ASCII_MATLAB);

private:
    void destroy()
    {
        if (A_ != nullptr)
        {
            MatDestroy(&A_);
        }
        A_ = nullptr;
    }

    /// PETSc matrix
    Mat A_ = nullptr;

    /// Number of the global rows
    PetscInt nrows_;

    /// Number of the global columns
    PetscInt ncols_;

    /// Number of the local rows
    PetscInt n_loc_rows_;

    /// Number of the local columns
    PetscInt n_loc_cols_;

    /// Starting index in a rank
    PetscInt start_rank_;

    /// Ending index in a rank
    PetscInt end_rank_;

    /*!
      \brief Create the matrix, configure memory allocation and set the
      related member data.
      \param d_nz Number of nonzeros per row in the diagonal portion of
                  local submatrix (same value is used for all local rows),
      \param o_nz Number of nonzeros per row in the off-diagonal portion of
                  local submatrix (same value is used for all local rows)
    */
    void create(const PetscInt d_nz, const PetscInt o_nz);

    friend bool finalizeMatrixAssembly(PETScMatrix& mat,
                                       const MatAssemblyType asm_type);
};

template <class T_DENSE_MATRIX>
void PETScMatrix::add(std::vector<PetscInt> const& row_pos,
                      std::vector<PetscInt> const& col_pos,
                      const T_DENSE_MATRIX& sub_mat)
{
    const PetscInt nrows = static_cast<PetscInt>(row_pos.size());
    const PetscInt ncols = static_cast<PetscInt>(col_pos.size());

    MatSetValues(A_, nrows, &row_pos[0], ncols, &col_pos[0], &sub_mat(0, 0),
                 ADD_VALUES);
};

/*!
    \brief          General interface for the matrix assembly.
    \param mat      The matrix to be finalized.
    \param asm_type Assembly type, either MAT_FLUSH_ASSEMBLY
                     or MAT_FINAL_ASSEMBLY.
*/
bool finalizeMatrixAssembly(
    PETScMatrix& mat, const MatAssemblyType asm_type = MAT_FINAL_ASSEMBLY);

}  // namespace MathLib
