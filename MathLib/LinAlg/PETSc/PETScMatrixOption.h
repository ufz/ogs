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

#include "petscmat.h"
#include "petscksp.h"

namespace MathLib
{
/*!
   \struct  PETScMatrixOption
   \brief This a struct data containing the configuration information to create a PETSc type matrix
*/
struct PETScMatrixOption
{
    PETScMatrixOption() :  _is_symmetric(false), _is_size_local_rows(false),
        _local_cols(PETSC_DECIDE), _d_nz(10), _o_nz(10)
    { }

    /// Flag for symmetric or unsymmetric matrix. The default is false.
    bool _is_symmetric;

    /*!
     \brief Flag for the type of the first argument of
            the constructor of class PETScMatrix, size.
            true:  the size is the number of local rows,
            false: the size is the number of global rows.
           The default is false.
     */
    bool _is_size_local_rows;

    /// Number of local columns. The default is PETSC_DECIDE.
    PetscInt _local_cols;

    /*!
     \brief Nnumber of nonzeros per row in DIAGONAL portion of local submatrix
           (same value is used for all local rows), the default is 10
    */
    PetscInt _d_nz;

    /*!
     \brief Number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix
            (same value is used for all local rows), the default is 10
    */
    PetscInt _o_nz;
};

} // end namespace
#endif

