/*!
   \file  PETScVector.cpp
   \brief Definition of member functions of class PETScVector,
          which provides an interface to PETSc vector routines.

     Note: the return message of PETSc routines is ommited in
           the source code. If it is really needed, it can be activated by
           adding a PetscErrorCode type variable before each PETSc fucntion

   \author Wenqing Wang
   \version
   \date Nov 2011 - Sep 2013

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScVector.h"

namespace MathLib
{
PETScVector::PETScVector(const PETScVector &existing_vec)
{
    _size = existing_vec._size;
    VecDuplicate(existing_vec._v, &_v);

    VecGetOwnershipRange(_v, &_start_rank,&_end_rank);
    VecGetLocalSize(_v, &_size_loc);
}

//-----------------------------------------------------------------
void PETScVector::create(PetscInt size, const PetscInt loc_size)
{
    if(loc_size == PETSC_DECIDE)
    {
        VecCreate(PETSC_COMM_WORLD, &_v);
        VecSetSizes(_v, loc_size, size);
    }
    else
    {
        // Fix size partitioning
        // the size can be associated to specific memory allocation of a matrix
        VecCreateMPI(PETSC_COMM_WORLD, loc_size, size, &_v);
    }
    VecSetFromOptions(_v);
    // VecSetUp(_v); // for ver.>3.3
    VecGetOwnershipRange(_v, &_start_rank, &_end_rank);

    VecGetLocalSize(_v, &_size_loc);
}

void PETScVector::gatherLocalVectors( PetscScalar local_array[],
                                      PetscScalar global_array[])
{
    // Collect vectors from processors.
    int size_rank;
    MPI_Comm_size(PETSC_COMM_WORLD, &size_rank);

    // number of elements to be sent for each rank
    std::vector<PetscInt>  i_cnt(size_rank);
    // offset in the receive vector of the data from each rank
    std::vector<PetscInt>  i_disp(size_rank);

    MPI_Allgather(&_size_loc, 1, MPI_INT, &i_cnt[0], 1, MPI_INT, PETSC_COMM_WORLD);

    // colloect local array
    PetscInt offset = 0;
    for(PetscInt i=0; i<size_rank; i++)
    {
        i_disp[i] = offset;
        offset += i_cnt[i];
    }

    MPI_Allgatherv(local_array, _size_loc, MPI_DOUBLE,
                   global_array, &i_cnt[0], &i_disp[0], MPI_DOUBLE, PETSC_COMM_WORLD);

}

void PETScVector::getGlobalVector(PetscScalar u[])
{

#ifdef TEST_MEM_PETSC
    PetscLogDouble mem1, mem2;
    PetscMemoryGetCurrentUsage(&mem1);
#endif

    PetscScalar *xp = nullptr;
    VecGetArray(_v, &xp);

    gatherLocalVectors(xp, u);

    //This following line may be needed late on
    //  for a communication load balance:
    //MPI_Barrier(PETSC_COMM_WORLD);

    VecRestoreArray(_v, &xp);

    //TEST
#ifdef TEST_MEM_PETSC
    PetscMemoryGetCurrentUsage(&mem2);
    PetscPrintf(PETSC_COMM_WORLD, "### Memory usage by Updating. Before :%f After:%f Increase:%d\n", mem1, mem2, (int)(mem2 - mem1));
#endif
}

PetscScalar PETScVector::getNorm(MathLib::VecNormType nmtype) const
{
    NormType petsc_norm = NORM_1;
    switch(nmtype)
    {
        case MathLib::VecNormType::NORM1:
            petsc_norm = NORM_1;
            break;
        case MathLib::VecNormType::NORM2:
            petsc_norm = NORM_2;
            break;
        case MathLib::VecNormType::INFINITY_N:
            petsc_norm = NORM_INFINITY;
            break;
        default:
            break;
    }

    PetscScalar norm = 0.;
    VecNorm(_v, petsc_norm, &norm);
    return norm;
}

void PETScVector::viewer(const std::string &file_name, const PetscViewerFormat vw_format)
{
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
    PetscViewerPushFormat(viewer, vw_format);

    finalizeAssembly();

    PetscObjectSetName((PetscObject)_v, file_name.c_str());
    VecView(_v, viewer);

#define  nEXIT_TEST
#ifdef EXIT_TEST
    VecDestroy(&_v);
    PetscFinalize();
    exit(0);
#endif

}

} //end of namespace

