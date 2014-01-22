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

#include "InfoMPI.h"

namespace MathLib
{

PETScVector::PETScVector(const PETScVector &existing_vec)
{
    _size = existing_vec._size;
    VecDuplicate(existing_vec._v, &_v);

    VecGetOwnershipRange(_v, &_start_rank,&_end_rank);
    VecGetLocalSize(_v, &_size_loc);

    // If values of the vector are copied too:
    //VecCopy(existing_vec._v, _v);

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
    VecGetOwnershipRange(_v, &_start_rank, &_end_rank);

    VecGetLocalSize(_v, &_size_loc);
}

void PETScVector::collectLocalVectors(  PetscScalar u_local_filled[],
                                        PetscScalar u_local_received[],
                                        PetscScalar u_global[])
{
    // Collect solution from processes.
    PetscInt receive_size_loc;
    PetscInt receive_start_rank;
    MPI_Status status;
    const int tag = 89999; // sending, receiving tag
    const int _size_rank = BaseLib::InfoMPI::getSize();
    const int _rank = BaseLib::InfoMPI::getRank();
    for(int i=0; i<_size_rank; i++)
    {
        if(i != _rank)
        {
            // get the local size belong to other ranks
            MPI_Sendrecv( &_size_loc, 1, MPI_INT, i,tag,
                          &receive_size_loc, 1, MPI_INT, i, tag, PETSC_COMM_WORLD, &status);
            // get the start point of the local vectors belong to other ranks
            MPI_Sendrecv( &_start_rank, 1, MPI_INT, i, tag,
                          &receive_start_rank, 1, MPI_INT, i, tag, PETSC_COMM_WORLD, &status);
            // get the other local vector
            MPI_Sendrecv( u_local_filled, _size_loc, MPI_DOUBLE, i,tag,
                          u_local_received, receive_size_loc, MPI_DOUBLE, i, tag, PETSC_COMM_WORLD, &status);
            for(int j=0; j<receive_size_loc; j++)
                u_global[receive_start_rank+j] = u_local_received[j];
        }
    }
}

void PETScVector::getGlobalEntries(PetscScalar u0[], PetscScalar u1[])
{

#ifdef TEST_MEM_PETSC
    PetscLogDouble mem1, mem2;
    PetscMemoryGetCurrentUsage(&mem1);
#endif

    PetscScalar *xp = nullptr;
    VecGetArray(_v, &xp);
    std::copy(xp, xp + _size_loc, u1);
    // Alternative for debugging:
    //for(int i=0; i<_size_loc; i++)
    // u1[i] = xp[i];

    PetscScalar *global_buff = new PetscScalar[_size];

    std::copy_n(u1, _size_loc, global_buff+_start_rank);
    // Alternative for debugging:
    //for(int j=0; j<_size_loc; j++)
    // global_buff[_start_rank+j] = u1[j];

    collectLocalVectors(u1, u0,  global_buff);

    //MPI_Barrier(PETSC_COMM_WORLD);
    // Copy the collected solution to the array for the previous solution
    for(int i=0; i<_size; i++)
    {
        u1[i] = global_buff[i];
        u0[i] = global_buff[i];
    }

    VecRestoreArray(_v, &xp);

    delete [] global_buff;

    //TEST
#ifdef TEST_MEM_PETSC
    PetscMemoryGetCurrentUsage(&mem2);
    PetscPrintf(PETSC_COMM_WORLD, "### Memory usage by Updating. Before :%f After:%f Increase:%d\n", mem1, mem2, (int)(mem2 - mem1));
#endif
}

PetscReal PETScVector::getNorm(OGS_NormType  nmtype) const
{
    NormType petsc_norm = NORM_1;
    switch(nmtype)
    {
        case OGS_NormType::SUM_ABS_ENTRIES:
            petsc_norm = NORM_1;
            break;
        case OGS_NormType::EUCLIDEAN:
            petsc_norm = NORM_2;
            break;
        case OGS_NormType::MAX_ABS_ENTRY:
            petsc_norm = NORM_INFINITY;
            break;
        default:
            break;
    }

    PetscReal norm = 0.;
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

