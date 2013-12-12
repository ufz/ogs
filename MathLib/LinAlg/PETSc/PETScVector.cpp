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

#include<iostream>
#include "BaseLib/MPI/InforMPI.h"

namespace MathLib
{

PETScVector::PETScVector(const PetscInt size)
{
   _size = size;
   create(_size);

   _size_loc = PETSC_DECIDE;
}

PETScVector::PETScVector(const PETScVector &existing_vec)
{
   _size = existing_vec._size;
   VecDuplicate(existing_vec._v, &_v);

   _size_loc = existing_vec._size_loc;

   VecGetOwnershipRange(_v, &_start_rank,&_end_rank);

   // If values of the vector are copied too:
   //VecCopy(existing_vec._v, _v);

}

PETScVector::~PETScVector()
{
   VecDestroy(&_v);
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
   VecGetOwnershipRange(_v, &_start_rank,&_end_rank);
}

void PETScVector::finalizeAssembly()
{
   VecAssemblyBegin(_v);
   VecAssemblyEnd(_v);
}

void PETScVector::getGlobalEntries(PetscScalar u0[], PetscScalar u1[])
{
#ifdef TEST_MEM_PETSC
   PetscLogDouble mem1, mem2;
   PetscMemoryGetCurrentUsage(&mem1);
#endif

   int i, j;
   PetscScalar *xp;

   int receivecount;
   PetscInt low,high,otherlow;
   MPI_Status status;
   PetscInt count;
   int tag = 89999; // sending, reveiving tag

   const int _size_rank = BaseLib::InforMPI::getSize();
   const int _rank = BaseLib::InforMPI::getRank();

   VecGetOwnershipRange(_v, &low, &high);
   VecGetLocalSize(_v, &count);

   VecGetArray(_v, &xp);
   for(i=0; i<count; i++)
      u1[i] = xp[i];

   PetscScalar *global_buff = new PetscScalar[_size];

   // Collect solution from processes.
   for(j=0; j<count; j++)
      global_buff[low+j] = u1[j];

   for(i=0; i<_size_rank; i++)
   {
      if(i != _rank)
      {
         MPI_Sendrecv( &count, 1, MPI_INT, i,tag,
                       &receivecount,1,MPI_INT,i,tag, PETSC_COMM_WORLD ,&status);
         MPI_Sendrecv( &low, 1, MPI_INT, i,tag,
                       &otherlow,1,MPI_INT,i,tag,PETSC_COMM_WORLD,&status );
         MPI_Sendrecv( u1, count, MPI_DOUBLE, i,tag,
                       u0,receivecount,MPI_DOUBLE,i,tag, PETSC_COMM_WORLD,&status  );
         for(j=0; j<receivecount; j++)
            global_buff[otherlow+j] = u0[j];
      }
   }

   //MPI_Barrier(PETSC_COMM_WORLD);
   // Copy the collected solution to the array for the previous solution
   for(i=0; i<_size; i++)
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

PetscInt PETScVector::getLocalVector(PetscScalar loc_vec[]) const
{
   PetscInt count;
   VecGetLocalSize(_v, &count);

   VecGetArray(_v, &loc_vec);

   return count;
}

void  PETScVector::getEntries(PetscInt ni,const PetscInt ix[],
                              PetscScalar y[]) const
{
   VecGetValues(_v, ni, ix, y);
}

PetscReal PETScVector::getNorm(VectorNormType  nmtype) const
{
   NormType nm_t[] = {NORM_1, NORM_2, NORM_INFINITY};

   PetscReal norm = 0.;
   VecNorm(_v, nm_t[nmtype], &norm);
   return norm;
}

void  PETScVector::restoreLocalVector(PetscScalar loc_vec[])
{
   VecRestoreArray(_v, &loc_vec);
}

void PETScVector::set(const int i, const PetscScalar value )
{

   VecSetValues(_v,1,&i,&value,INSERT_VALUES);
}

void PETScVector::add(const int i, const PetscScalar value)
{
   VecSetValue(_v, i, value,  ADD_VALUES);
}

void  PETScVector::setValues( PetscInt ni, const PetscInt ix[],
                              const PetscScalar y[])
{
   VecSetValues(_v, ni, ix, y, INSERT_VALUES);
}

void  PETScVector::addValues( PetscInt ni, const PetscInt ix[],
                              const PetscScalar y[])
{
   VecSetValues(_v, ni, ix, y, ADD_VALUES);
}

void PETScVector::setZero( )
{
   VecSet(_v, 0.0);
}

double  PETScVector::get(const  PetscInt idx) const
{
   double x[1];
   PetscInt idxs[1];
   idxs[0] = idx;

   VecGetValues(_v, 1, idxs, x);
   return x[0];
}

// Overloaded operator: initialize  the vector with a constant value
void PETScVector::operator= (const PetscScalar val)
{
   VecSet(_v, val);
}

//Overloaded operator: assignment
void PETScVector::operator= (PETScVector &v_in)
{
   VecCopy(_v, v_in._v);
}

//Overloaded operator:add
void PETScVector::operator+= (const PETScVector& v_in)
{
   VecAXPY(_v, 1.0, v_in._v);
}

// Overloaded operator: subtract
void PETScVector::operator-= (const PETScVector& v_in)
{
   VecAXPY(_v, -1.0, v_in._v);
}

void PETScVector::Viewer(const std::string &file_name)
{
   PetscViewer viewer;
   const std::string fname = file_name + "_petsc_global_vector_entries.txt";
   PetscViewerASCIIOpen(PETSC_COMM_WORLD, fname.c_str(), &viewer);
   // PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
   PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_VTK);

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

