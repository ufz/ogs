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

PETScMatrix::PETScMatrix (const PetscInt nrows, const PETScMatrixOption &mat_opt)
    :_nrows(nrows), _ncols(nrows), _n_loc_rows(PETSC_DECIDE),
     _n_loc_cols(mat_opt.n_local_cols)
{
    if(!mat_opt.is_global_size)
    {
        _nrows = PETSC_DECIDE;
        _ncols = PETSC_DECIDE;
        _n_loc_rows = nrows;

        // Make the matrix be square.
        MPI_Allreduce(&_n_loc_rows, &_nrows, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
        _ncols = _nrows;
    }

    create(mat_opt.d_nz, mat_opt.o_nz);
}

PETScMatrix::PETScMatrix (const PetscInt nrows, const PetscInt ncols, const PETScMatrixOption &mat_opt)
    :_nrows(nrows), _ncols(ncols),  _n_loc_rows(PETSC_DECIDE),
     _n_loc_cols(mat_opt.n_local_cols)
{
    if(!mat_opt.is_global_size)
    {
        _nrows = PETSC_DECIDE;
        _ncols = PETSC_DECIDE;
        _n_loc_rows = nrows;
        _n_loc_cols = ncols;
    }

    create(mat_opt.d_nz, mat_opt.o_nz);
}

void PETScMatrix::setRowsColumnsZero(std::vector<PetscInt> const& row_pos)
{
    // Each rank (compute core) processes only the rows that belong to the rank itself.
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

#define nEXIT_TEST
#ifdef EXIT_TEST
    MatDestroy(&_A);
    PetscFinalize();
    exit(0);
#endif

}

void PETScMatrix::create(const PetscInt d_nz, const PetscInt o_nz)
{
    MatCreate(PETSC_COMM_WORLD, &_A);
    MatSetSizes(_A, _n_loc_rows, _n_loc_cols, _nrows, _ncols);

    MatSetFromOptions(_A);

    // for a dense matrix: MatSeqAIJSetPreallocation(_A, d_nz, PETSC_NULL);
    MatMPIAIJSetPreallocation(_A, d_nz, PETSC_NULL, o_nz, PETSC_NULL);

    MatGetOwnershipRange(_A, &_start_rank, &_end_rank);
    MatGetSize(_A, &_nrows,  &_ncols);
    MatGetLocalSize(_A, &_n_loc_rows, &_n_loc_cols);
}

bool finalizeMatrixAssembly(PETScMatrix &mat, const MatAssemblyType asm_type)
{
    mat.finalizeAssembly(asm_type);
    return true;
}

} //end of namespace

