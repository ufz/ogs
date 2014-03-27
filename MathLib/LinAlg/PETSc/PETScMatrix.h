/*!
   \file  PETScMatrix.h
   \brief Declaration of class PETScMatrix, which provides an interface to
          PETSc matrix routines.

   \author Wenqing Wang
   \date Nov 2013 - 2014

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef PETSCMATRIX_H_
#define PETSCMATRIX_H_

#include <string>
#include <vector>

#include "PETScMatrixOption.h"
#include "PETScVector.h"

typedef Mat PETSc_Mat;

namespace MathLib
{

/*!
   \brief Wrapper class for PETSc matrix routines for a square matrix for
          the system of linear equations.
*/
class PETScMatrix
{
    public:
        /*!
          \brief       Constructor for the matrix partitioning with default options:
                       The size of the glabal size, the numbers of local rows
                       and columns have the value of PETSC_DECIDE.
          \param size  The dimension of the matrix.
        */
        explicit PETScMatrix(const PetscInt size);

        /*!
          \brief        Constructor for the user determined partitioning with more options
          \param size   The number of rows of the matrix or the local matrix.
          \param mat_op The configuration information for creating a matrix.
        */
        PETScMatrix(const PetscInt size, const PETScMatrixOption &mat_op = PETScMatrixOption() );

        ~PETScMatrix()
        {
            MatDestroy(&_A);
        }

        /*!
           \brief          Perform MPI collection of assembled entries in buffer
           \param asm_type Assmebly type, either MAT_FLUSH_ASSEMBLY
                           or MAT_FINAL_ASSEMBLY
        */
        void finalizeAssembly(const MatAssemblyType asm_type = MAT_FINAL_ASSEMBLY)
        {
            MatAssemblyBegin(_A, asm_type);
            MatAssemblyEnd(_A, asm_type);
        }

        /// Get the number of rows.
        PetscInt getNRows() const
        {
            return _size;
        }

        /// Get the number of columns.
        PetscInt getNCols() const
        {
            return _size;
        }

        /// Get the start global index of the rows of the same rank.
        PetscInt getRangeBegin() const
        {
            return _start_rank;
        }

        /// Get the end global index of the rows in the same rank.
        PetscInt getRangeEnd() const
        {
            return _end_rank;
        }

        /// Get the number of local rows.
        PetscInt getNLocalRows() const
        {
            return _n_loc_rows;
        }

        /// Get the number of local columns.
        PetscInt getNLocalColumns() const
        {
            return _n_loc_cols;
        }

        /// Get matrix reference.
        PETSc_Mat &getRawMatrix()
        {
            return _A;
        }

        /// Set all entries to zero.
        void setZero()
        {
            MatZeroEntries(_A);
        }

        /*
           \brief Set the specified rows to zero except diagonal entries, i.e.
                  \f$A(k, j) = 0, j!=k, j=1,2,\cdots, n\f$, where \f$k \in \mbox{row\_pos}\f$
                  This fucntion must be called by all rank.
           \param row_pos The row indicies of the specified rows.
        */
        void setRowsColumnsZero(std::vector<PetscInt> const& row_pos);

        /*
           \brief       Perform operation \f$ y = A x \f$
           \param vec   The given vector, e.g. \f$ x \f$
           \param vec_r The result vector, e.g. \f$ y \f$
            Both of the two arguments must be created prior to be used.
        */
        void multi(const PETScVector &vec, PETScVector &vec_r)
        {
            MatMult(_A, vec.getData(), vec_r.getData() );
        }

        /*!
           \brief       Set a single entry with a value.
           \param i     The row index.
           \param j     The column index.
           \param value The entry value.
        */
        void set(const PetscInt i, const PetscInt j, const PetscScalar value)
        {
            MatSetValue(_A, i, j, value, INSERT_VALUES);
        }

        /*!
           \brief Add value to a single entry.
           \param i     The row index.
           \param j     The column index.
           \param value The entry value.
        */
        void add(const PetscInt i, const PetscInt j, const PetscScalar value)
        {
            MatSetValue(_A, i, j, value, ADD_VALUES);
        }

        /*!
          \brief         Add a submatrix to this.
          \param row_pos The row indices of the entries of the submatrix.
          \param col_pos The column indices of the entries of the submatrix.
          \param sub_mat A dense matrix to be added on.
        */
        template <class T_DENSE_MATRIX>
        void add(std::vector<PetscInt> const& row_pos,
                 std::vector<PetscInt> const& col_pos,
                 const T_DENSE_MATRIX &sub_mat );

        /*! View the global vector for test purpose. Do not use it for output a big vector.
            \param file_name  File name for output
            \vw_format        File format listed as:
             PETSC_VIEWER_DEFAULT            Default format
             PETSC_VIEWER_ASCII_MATLAB       MATLAB format
             PETSC_VIEWER_ASCII_DENSE        Print matrix as dense
             PETSC_VIEWER_ASCII_IMPL         Implementation-specific format
                                               (which is in many cases the same as the default)
             PETSC_VIEWER_ASCII_INFO         Basic information about object
             PETSC_VIEWER_ASCII_INFO_DETAIL  More detailed info about object
             PETSC_VIEWER_ASCII_COMMON       Identical output format for all objects of a particular type
             PETSC_VIEWER_ASCII_INDEX        (for vectors) Prints the vector element number next to
                                                each vector entry
             PETSC_VIEWER_ASCII_SYMMODU      Print parallel vectors without indicating the processor ranges
             PETSC_VIEWER_ASCII_VTK          Outputs the object to a VTK file
             PETSC_VIEWER_NATIVE             Store the object to the binary file in its native format
                                              (for example, dense matrices are stored as dense),
                                              DMDA vectors are dumped directly to the file instead of
                                              being first put in the natural ordering
             PETSC_VIEWER_DRAW_BASIC         Views the vector with a simple 1d plot
             PETSC_VIEWER_DRAW_LG            Views the vector with a line graph
             PETSC_VIEWER_DRAW_CONTOUR       Views the vector with a contour plot
        */
        void viewer(const std::string &file_name,
                    const PetscViewerFormat vw_format = PETSC_VIEWER_ASCII_MATLAB );

    private:
        /// PETSc matrix
        PETSc_Mat _A;

        /// Dimension of matrix
        PetscInt _size;
        /// Number of the local rows
        PetscInt _n_loc_rows;
        /// Number of the local columns
        PetscInt _n_loc_cols;

        /// Starting index in a rank
        PetscInt _start_rank;
        /// Ending index in a rank
        PetscInt _end_rank;

        /// Create the matrix
        void create();

        /*!
          \brief Config memory allocation and set the related member data.
          \param Number of nonzeros per row in the diagonal portion of local submatrix
                 (same value is used for all local rows),
          \param Number of nonzeros per row in the off-diagonal portion of local submatrix
                 (same value is used for all local rows)
        */
        void config(const PetscInt d_nz, const PetscInt o_nz);

        friend bool finalizeMatrixAssembly(PETScMatrix &mat, const MatAssemblyType asm_type);
};

/*!
    \brief          Add a dense sub-matrix to a PETSc matrix.
    \param row_pos  The global indices of the rows of the dense sub-matrix.
    \param col_pos  The global indices of the colums of the dense sub-matrix.
    \param sub_mat  A dense sub-matrix to be added.
*/
template<class T_DENSE_MATRIX>
void PETScMatrix::add(std::vector<PetscInt> const& row_pos,
                      std::vector<PetscInt> const& col_pos,
                      const T_DENSE_MATRIX &sub_mat)
{
    const PetscInt nrows = static_cast<PetscInt> (row_pos.size());
    const PetscInt ncols = static_cast<PetscInt> (col_pos.size());

    MatSetValues(_A, nrows, &row_pos[0], ncols, &col_pos[0], sub_mat.getEntries(), ADD_VALUES);
};

/*!
    \brief          General interface for the matrix assembly.
    \param mat      The matrix to be finalized.
    \param asm_type Assmebly type, either MAT_FLUSH_ASSEMBLY
                     or MAT_FINAL_ASSEMBLY.
*/
bool finalizeMatrixAssembly(PETScMatrix &mat, const MatAssemblyType asm_type = MAT_FINAL_ASSEMBLY);

} // end namespace
#endif

