/*!
   \file  PETScMatrixOption.h
   \brief Define data for the configuration of PETSc matrix and linear solver.

   \author Wenqing Wang
   \date 02-2014

   \copyright
    Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#pragma once

#include <petscmat.h>

namespace MathLib
{
/*!
   \brief This a struct data containing the configuration information to create
   a PETSc type matrix
*/
struct PETScMatrixOption
{
    PETScMatrixOption()
        : is_global_size(true),
          n_local_cols(PETSC_DECIDE),
          d_nz(PETSC_DECIDE),
          o_nz(PETSC_DECIDE)
    {
    }

    /*!
     \brief Flag for the type of size, which is one of arguments of
            the constructor of class PETScMatrix
              true:  the size is the number of local rows,
              false: the size is the number of global rows.
            The default is false.
     */
    bool is_global_size;

    /// Number of local columns. The default is PETSC_DECIDE.
    PetscInt n_local_cols;

    /*!
     \brief Number of nonzeros per row in the diagonal portion of local
     submatrix
           (same value is used for all local rows), the default is PETSC_DECIDE
    */
    PetscInt d_nz;

    /*!
     \brief Number of nonzeros per row in the off-diagonal portion of local
     submatrix
            (same value is used for all local rows), the default is PETSC_DECIDE
    */
    PetscInt o_nz;
};

}  // end namespace
