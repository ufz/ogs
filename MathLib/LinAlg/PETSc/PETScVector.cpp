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

#include "array"

#include "InfoMPI.h"

namespace MathLib
{
PETScVector::PETScVector(const PETScVector &existing_vec)
{
    Duplicate(existing_vec);
    VecCopy(existing_vec._v, _v);
}

void PETScVector::Duplicate(const PETScVector &existing_vec)
{
    _size = existing_vec._size;
    VecDuplicate(existing_vec._v, &_v);

    VecGetOwnershipRange(_v, &_start_rank,&_end_rank);
    VecGetLocalSize(_v, &_size_loc);
}

//-----------------------------------------------------------------
void PETScVector::create(PetscInt vec_size)
{
    VecCreate(PETSC_COMM_WORLD, &_v);
    // The following two lines are used to test a fix size partition
    // VecCreateMPI(PETSC_COMM_WORLD,m_size_loc, m, &v);
    // VecSetSizes(v, m_size_loc, m);
    VecSetSizes(_v, PETSC_DECIDE, vec_size);
    VecSetFromOptions(_v);
    // VecSetUp(_v); // for ver.>3.3
    VecGetOwnershipRange(_v, &_start_rank, &_end_rank);

    VecGetLocalSize(_v, &_size_loc);
}

void PETScVector::collectLocalVectors( PetscScalar local_array[],
                                       PetscScalar global_array[])
{
    // Collect solution from processes.
    const int size_rank = BaseLib::InfoMPI::getSize();

    // number of elements to be sent for each rank
    int *i_cnt = new int[size_rank];
    // offset in the receive vector of the data from each rank
    int *i_disp = new int[size_rank];

    MPI_Allgather(&_size_loc, 1, MPI_INT, i_cnt, 1, MPI_INT, PETSC_COMM_WORLD);

    // colloect local array
    int offset = 0;
    for(int i=0; i<size_rank; i++)
    {
        i_disp[i] = offset;
        offset += i_cnt[i];
    }

    MPI_Allgatherv(local_array, _size_loc, MPI_DOUBLE,
                   global_array, i_cnt, i_disp, MPI_DOUBLE, PETSC_COMM_WORLD);

    delete [] i_cnt;
    delete [] i_disp;
}

void PETScVector::getGlobalEntries(PetscScalar u[])
{

#ifdef TEST_MEM_PETSC
    PetscLogDouble mem1, mem2;
    PetscMemoryGetCurrentUsage(&mem1);
#endif

    PetscScalar *xp = nullptr;
    VecGetArray(_v, &xp);

    collectLocalVectors(xp, u);

    //MPI_Barrier(PETSC_COMM_WORLD);

    VecRestoreArray(_v, &xp);

    //TEST
#ifdef TEST_MEM_PETSC
    PetscMemoryGetCurrentUsage(&mem2);
    PetscPrintf(PETSC_COMM_WORLD, "### Memory usage by Updating. Before :%f After:%f Increase:%d\n", mem1, mem2, (int)(mem2 - mem1));
#endif
}

PetscScalar PETScVector::getNorm(MathLib::VecNormType  nmtype) const
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

PetscScalar  PETScVector::get(const  PetscInt idx) const
{
    double x;
    VecGetValues(_v, 1, &idx, &x);
    return x;
}

void PETScVector::viewer(const std::string &file_name, const PetscViewerFormat vw_format)
{
    PetscViewer viewer;
    const std::string fname = file_name + "_petsc_global_vector_entries.txt";
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, fname.c_str(), &viewer);
    PetscViewerPushFormat(viewer, vw_format);

    finalizeAssembly();

    PetscObjectSetName((PetscObject)_v,file_name.c_str());
    VecView(_v, viewer);

#define  nEXIT_TEST
#ifdef EXIT_TEST
    VecDestroy(&_v);
    PetscFinalize();
    exit(0);
#endif

}

} //end of namespace

