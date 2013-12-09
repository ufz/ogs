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


namespace MathLib
{


PETScVector :: PETScVector ()
{
   _size_loc = PETSC_DECIDE;
}

PETScVector:: PETScVector(const PetscInt size)
{
   _size = size;
   create(_size);

   _size_loc = PETSC_DECIDE;
}

PETScVector:: PETScVector(const PETScVector &existing_vec)
{

   _size = existing_vec._size;
   VecDuplicate(existing_vec.v, &v);

   _size_loc = existing_vec._size_loc;

   VecGetOwnershipRange(v, &_start_rank,&_end_rank);

   // If values of the vector are copied too:
   //VecCopy(existing_vec.v, v);

}

PETScVector:: ~PETScVector()
{
   VecDestroy(&v);
}

void PETScVector::init(const PetscInt vec_size)
{
   _size = vec_size;
   create(_size);
}

//-----------------------------------------------------------------
void  PETScVector::create(PetscInt vec_size)
{
   VecCreate(PETSC_COMM_WORLD, &v);
   // The following two lines are used to test a fix size partition
   // VecCreateMPI(PETSC_COMM_WORLD,m_size_loc, m, &v);
   // VecSetSizes(v, m_size_loc, m);
   VecSetSizes(v, PETSC_DECIDE, vec_size);
   VecSetFromOptions(v);
   VecGetOwnershipRange(v, &_start_rank,&_end_rank);
}


void PETScVector::finalizeAssembly()
{
   VecAssemblyBegin(v);
   VecAssemblyEnd(v);
}

void PETScVector::getGlobalEntries(PetscScalar u0[], PetscScalar u1[])
{

#ifdef TEST_MEM_PETSC
   //TEST
   PetscLogDouble mem1, mem2;
   PetscMemoryGetCurrentUsage(&mem1);
#endif

   int i, j;
   PetscScalar *xp;

   int receivecount;
   PetscInt low,high,otherlow;
   MPI_Status status;
   PetscInt count;
   int tag = 9999;
   VecGetOwnershipRange(v, &low, &high);
   VecGetLocalSize(v, &count);

   VecGetArray(v, &xp);
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

   VecRestoreArray(v, &xp);

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
   VecGetLocalSize(v, &count);

   VecGetArray(v, &loc_vec);

   return count;
}

void  PETScVector::getEntries(PetscInt ni,const PetscInt ix[],
                              PetscScalar y[]) const
{
   VecGetValues(v, ni, ix, y);
}

PetscReal PETScVector::getNorm(NormType  nmtype) const
{
   PetscReal norm = 0.;
   VecNorm(v, nmtype, &norm);
   return norm;
}

void  PETScVector::restoreLocalVector(PetscScalar loc_vec[])
{
   VecRestoreArray(v, &loc_vec);
}

void PETScVector::set(const int i, const PetscScalar value )
{

   VecSetValues(v,1,&i,&value,INSERT_VALUES);
}

void  PETScVector::setValues( PetscInt ni, const PetscInt ix[],
                              const PetscScalar y[],InsertMode mode)
{
   VecSetValues(v, ni, ix, y, mode);
}

void PETScVector::add(const int i, const PetscScalar value,InsertMode mode )
{

   VecSetValue(v, i, value, mode);
}

void PETScVector::setZero( )
{

   VecSet(v, 0.0);
}

double  PETScVector::get(const  PetscInt idx) const
{
   double x[1];
   PetscInt idxs[1];
   idxs[0] = idx;

   VecGetValues(v, 1, idxs, x);
   return x[0];
}


// Overloaded operator: initialize  the vector with a constant value
void PETScVector::operator= (const PetscScalar val)
{

   VecSet(v, val);
}

//Overloaded operator: assignment
PETScVector& PETScVector::operator= (PETScVector &v_in)
{
   VecCopy(v, v_in.v);
   return *this;
}

//Overloaded operator:add
void PETScVector::operator+= (const PETScVector& v_in)
{
   VecAXPY(v, 1.0, v_in.v);
}

// Overloaded operator: subtract
void PETScVector::operator-= (const PETScVector& v_in)
{
   VecAXPY(v, -1.0, v_in.v);
}

void PETScVector::Viewer(const std::string &file_name)
{
   PetscViewer viewer;
   const std::string fname = file_name + "_petsc_global_vector_entries.txt";
   PetscViewerASCIIOpen(PETSC_COMM_WORLD, fname.c_str(), &viewer);
   PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);

   finalizeAssembly();

   // PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_VTK);
   PetscObjectSetName((PetscObject)v,file_name.c_str());
   VecView(v, viewer);

#define  nEXIT_TEST
#ifdef EXIT_TEST
   VecDestroy(&v);
   PetscFinalize();
   exit(0);
#endif

}

} //end of namespace

