/*!
   \file  PETScMatrix.cpp
   \brief Definition of member functions of class PETScMatrix, which provides an interface to
          PETSc matrix routines.

   \author Wenqing Wang
   \date Nov 2011 - Sep 2013

   \copyright
   Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScMatrix.h"

namespace MathLib
{

PETScMatrix::PETScMatrix (const PetscInt size, const PETScMatrixOption mat_opt)
{
    _size = size;
    _loc_rows = PETSC_DECIDE;
    _loc_cols = mat_opt._local_cols;

    if(mat_opt._is_size_local_rows)
    {
        _size = PETSC_DECIDE;
        _loc_rows = size;
    }

    MatCreate(PETSC_COMM_WORLD, &_A);
    MatSetSizes(_A, _loc_rows, _loc_cols, _size, _size);

    MatSetFromOptions(_A);

    // for a dense matrix: MatSeqAIJSetPreallocation(_A, d_nz, PETSC_NULL);
    MatMPIAIJSetPreallocation(_A, mat_opt._d_nz, PETSC_NULL, mat_opt._o_nz, PETSC_NULL);

    MatGetOwnershipRange(_A, &_start_rank, &_end_rank);
    MatGetSize(_A, &_size,  PETSC_NULL);
    MatGetLocalSize(_A, &_loc_rows, &_loc_cols);
}

//-----------------------------------------------------------------
void PETScMatrix::setRowsColumnsZero(std::vector<PetscInt> const& row_pos)
{
    // Each process indicates only rows it owns that are to be zeroed.
    // If it is called, it must be called by all ranks of cores.

    const PetscScalar one = 1.0;
    const PetscInt nrows = static_cast<PetscInt> (row_pos.size());

    if(nrows>0)
        MatZeroRows(_A, nrows, &row_pos[0], one, PETSC_NULL, PETSC_NULL);
    else
        MatZeroRows(_A, 0, PETSC_NULL, one, PETSC_NULL, PETSC_NULL);
}

void PETScMatrix::viewer(const std::string &file_name, const PetscViewerFormat vw_format)
{
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
    PetscViewerPushFormat(viewer, vw_format);

    finalizeAssembly();

    PetscObjectSetName((PetscObject)_A,"Stiffness_matrix");
    MatView(_A,viewer);

#define  nEXIT_TEST
#ifdef EXIT_TEST
    MatDestroy(&_A);
    PetscFinalize();
    exit(0);
#endif

}

bool finalizeMatrixAssembly(PETScMatrix &mat, const MatAssemblyType asm_type)
{
    mat.finalizeAssembly(asm_type);
    return true;
}

} //end of namespace

