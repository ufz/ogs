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
   \class PETScVector

   \brief Wrapper class for PETSc matrix routines
*/
class PETScMatrix
{
    public:
        /*!
          \brief Constructor for the PETSc determined partitioning.
          \param size   The dimension of the matrix.
          \param mat_op The configuration information for creating a matrix
        */
        PETScMatrix(const PetscInt size, const PETScMatrixOption mat_op = PETScMatrixOption() );

        ~PETScMatrix()
        {
            MatDestroy(&_A);
        }

        /*!
           \brief Perform MPI collection of assembled entries in buffer
           \param asm_type Assmebly type, either MAT_FLUSH_ASSEMBLY
                           or MAT_FINAL_ASSEMBLY
        */
        void finalizeAssembly(const MatAssemblyType asm_type = MAT_FINAL_ASSEMBLY)
        {
            MatAssemblyBegin(_A, asm_type);
            MatAssemblyEnd(_A, asm_type);
        }

        /// Set the matrix being symmetric
        void enableSymmetric()
        {
            finalizeAssembly();
            MatSetOption(_A, MAT_SYMMETRIC, PETSC_TRUE);
        }

        /// Get the dimension
        PetscInt size() const
        {
            return _size;
        }

        /// Get the start global index of the rows of the same rank
        PetscInt getRangeBegin() const
        {
            return _start_rank;
        }

        /// Get the last global index of the row in the same rank
        PetscInt getRangeEnd() const
        {
            return _end_rank;
        }

        /// Get the number of local rows
        PetscInt getLocalRows() const
        {
            return _loc_rows;
        }

        /// Get the number of local rows
        PetscInt getLocalColumns() const
        {
            return _loc_cols;
        }

        /// Is the matrix symmetric
        bool isSymmetric() const
        {
            return _is_symmetric;
        }

        /// Get matrix reference
        PETSc_Mat &getData()
        {
            return _A;
        }

        /// Set all entries to zero
        void setZero()
        {
            MatZeroEntries(_A);
        }

        /*
           \brief Set the specified rows and columns to zero except off-diagonal entries
           \param row_pos The row indicies of the specified rows.
        */
        void setRowsColumnsZero(std::vector<PetscInt> const& row_pos);

        /*
           \brief Perform operation \f$ y = A x \f$
           \param vec   The given vector, e.g. \f$ x \f$
           \param vec_r The result vector, e.g. \f$ y \f$
            Both of the two arguments must be created prior to be used.
        */
        void multVector(PETScVector &vec, PETScVector &vec_r)
        {
            MatMult(_A, vec.getData(), vec_r.getData() );
        }

        /*!
           \brief Insert a single entry with value.
           \param i     The row index
           \param j     The column index
           \param value The entry value
        */
        void set(const PetscInt i, const PetscInt j, const PetscScalar value)
        {
            MatSetValue(_A, i, j, value, INSERT_VALUES);
        }

        /*!
           \brief Add value to a single entry.
           \param i     The row index
           \param j     The column index
           \param value The entry value
        */
        void add(const PetscInt i, const PetscInt j, const PetscScalar value)
        {
            MatSetValue(_A, i, j, value, ADD_VALUES);
        }

        /*!
          \brief Add a submatrix to this

          \param row_pos The global row indicies of the entries of the submatrix
          \param col_pos The global column indicies of the entries of the submatrix
          \param sub_mat A dense matrix to be added on
        */
        template <class T_DENSE_MATRIX>
        void add(std::vector<PetscInt> const& row_pos,
                 std::vector<PetscInt> const& col_pos,
                 const T_DENSE_MATRIX &sub_mat );

        /*! View the global vector for test purpose. Do not use it for output a big vector.
            \param file_name  File name for output
            \vw_format        File format listed as:
        PETSC_VIEWER_DEFAULT 	- default format
        PETSC_VIEWER_ASCII_MATLAB 	- MATLAB format
        PETSC_VIEWER_ASCII_DENSE 	- print matrix as dense
        PETSC_VIEWER_ASCII_IMPL 	- implementation-specific format (which is in many cases the same as the default)
        PETSC_VIEWER_ASCII_INFO 	- basic information about object
        PETSC_VIEWER_ASCII_INFO_DETAIL 	- more detailed info about object
        PETSC_VIEWER_ASCII_COMMON 	- identical output format for all objects of a particular type
        PETSC_VIEWER_ASCII_INDEX 	- (for vectors) prints the vector element number next to each vector entry
        PETSC_VIEWER_ASCII_SYMMODU 	- print parallel vectors without indicating the processor ranges
        PETSC_VIEWER_ASCII_VTK 	- outputs the object to a VTK file
        PETSC_VIEWER_NATIVE 	- store the object to the binary file in its native format (for example, dense matrices are stored as dense), DMDA vectors are dumped directly to the file instead of being first put in the natural ordering
        PETSC_VIEWER_DRAW_BASIC 	- views the vector with a simple 1d plot
        PETSC_VIEWER_DRAW_LG 	- views the vector with a line graph
        PETSC_VIEWER_DRAW_CONTOUR 	- views the vector with a contour plot

        */
        void viewer(const std::string &file_name,
                    const PetscViewerFormat vw_format = PETSC_VIEWER_ASCII_MATLAB );

    private:
        /// PETSc matrix
        PETSc_Mat _A;

        /// Dimension
        PetscInt _size;
        /// Number of the local rows
        PetscInt _loc_rows;
        /// Number of the local columns
        PetscInt _loc_cols;

        /// Starting index in a rank
        PetscInt _start_rank;
        /// Ending index in a rank
        PetscInt _end_rank;

        /// Flag for symmetric or unsymmetric matrix. The default is false.
        bool _is_symmetric;

        friend bool finalizeMatrixAssembly(PETScMatrix &mat, const MatAssemblyType asm_type);
};

template<class T_DENSE_MATRIX>
void PETScMatrix::add(std::vector<PetscInt> const& row_pos,
                      std::vector<PetscInt> const& col_pos,
                      const T_DENSE_MATRIX &sub_mat)
{
    const PetscInt nrows = static_cast<PetscInt> (row_pos.size());
    const PetscInt ncols = static_cast<PetscInt> (col_pos.size());

    MatSetValues(_A, nrows, &row_pos[0], ncols, &col_pos[0], sub_mat.getData(), ADD_VALUES);
};

/*!
    \brief General interface for the matrix assembly
    \param mat      The matrix to be finalized
    \param asm_type Assmebly type, either MAT_FLUSH_ASSEMBLY
                     or MAT_FINAL_ASSEMBLY
*/
bool finalizeMatrixAssembly(PETScMatrix &mat, const MatAssemblyType asm_type = MAT_FINAL_ASSEMBLY);

} // end namespace
#endif

