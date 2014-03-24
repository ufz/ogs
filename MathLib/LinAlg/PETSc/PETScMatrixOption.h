/*!
   \file  PETScMatrixOption.h
   \brief Define data for the configuration of PETSc matrix and linear solver.

   \author Wenqing Wang
   \date 02-2014

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef PETSCMATRIXOPTION_H_
#define PETSCMATRIXOPTION_H_

#include <petscmat.h>

namespace MathLib
{
/*!
   \brief This a struct data containing the configuration information to create a PETSc type matrix
*/
struct PETScMatrixOption
{
    PETScMatrixOption() :  _is_global_size(true), _n_local_cols(PETSC_DECIDE),
        _d_nnz(PETSC_DECIDE), _o_nnz(PETSC_DECIDE)
    { }

    /*!
     \brief Flag for the type of size, which is one of arguments of
            the constructor of class PETScMatrix
              true:  the size is the number of local rows,
              false: the size is the number of global rows.
            The default is false.
     */
    bool _is_global_size;

    /// Number of local columns. The default is PETSC_DECIDE.
    PetscInt _n_local_cols;

    /*!
     \brief Number of nonzeros per row in DIAGONAL portion of local submatrix
           (same value is used for all local rows), the default is PETSC_DECIDE
    */
    PetscInt _d_nnz;

    /*!
     \brief Number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix
            (same value is used for all local rows), the default is PETSC_DECIDE
    */
    PetscInt _o_nnz;
};

} // end namespace
#endif

