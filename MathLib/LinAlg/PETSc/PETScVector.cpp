/*!
   \file  PETScVector.cpp
   \brief Definition of member functions of class PETScVector, which provides an interface
          to PETSc vector routines.

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
   m_size_loc = PETSC_DECIDE;
 }

PETScVector:: PETScVector(const PetscInt size)
{
   m_size = size;    
   Create(m_size); 

   m_size_loc = PETSC_DECIDE; 
}



PETScVector:: ~PETScVector()
{
  VecDestroy(&v);
}

void PETScVector::Init(const PetscInt size)
{
   m_size = size;    
   Create(m_size);   
}

//-----------------------------------------------------------------
void  PETScVector::Create(PetscInt m)
{
  //PetscErrorCode ierr;  // returned value from PETSc functions 
  VecCreate(PETSC_COMM_WORLD, &v);
  ////VecCreateMPI(PETSC_COMM_WORLD,m_size_loc, m, &v);
  //VecSetSizes(v, m_size_loc, m);
  VecSetSizes(v, PETSC_DECIDE, m);
  VecSetFromOptions(v);
  //VecGetOwnershipRange(v, &i_start,&i_end);
}

void PETScVector::CloneMe(PETScVector &new_vec)
{
    VecDuplicate(v, &new_vec.v);
}

void  PETScVector::getOwnerRange(int *start_r, int *end_r)
{
  *start_r = i_start;
  *end_r = i_end;
}


  void PETScVector::finalAssemble()
{
  VecAssemblyBegin(v);
  VecAssemblyEnd(v);
}


void PETScVector::getGlobalEntries(PetscScalar *u0, PetscScalar *u1)
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


  PetscScalar *global_buff = new double[m_size];


  // Collect solution from processes.
  for(j=0; j<count; j++)
    global_buff[low+j] = u1[j];
  for(i=0;i<mpi_size;i++) 
  {
     if(i != rank)
     {

       MPI_Sendrecv( &count, 1, MPI_INT, i,tag, 
                     &receivecount,1,MPI_INT,i,tag, PETSC_COMM_WORLD ,&status);
       MPI_Sendrecv( &low, 1, MPI_INT, i,tag,
                 &otherlow,1,MPI_INT,i,tag,PETSC_COMM_WORLD,&status );
       MPI_Sendrecv( u1, count, MPI_DOUBLE, i,tag,
                     u0,receivecount,MPI_DOUBLE,i,tag, PETSC_COMM_WORLD,&status  );
       for(j=0;j<receivecount;j++)
         global_buff[otherlow+j] = u0[j];
     }
  }


  //MPI_Barrier(PETSC_COMM_WORLD);
  // Copy the collected solution to the array for the previous solution
  for(i=0;i<m_size;i++)
  {
    u1[i] = global_buff[i];
    u0[i] = global_buff[i];
  }
 

  delete [] global_buff;

  VecRestoreArray(v, &xp);


  //TEST
#ifdef TEST_MEM_PETSC
   PetscMemoryGetCurrentUsage(&mem2);
   PetscPrintf(PETSC_COMM_WORLD, "### Memory usage by Updating. Before :%f After:%f Increase:%d\n", mem1, mem2, (int)(mem2 - mem1));
#endif

}



int PETScVector::getLocalVector(PetscScalar *loc_vec)
{
  PetscInt count;
  VecGetLocalSize(v, &count);

  VecGetArray(v, &loc_vec);

  return count;
}



/*!
  Get values of the specified elements from a global vector

  @param v_type - Indicator for vector: 0: x; 1: rhs
  @param ni 	- number of elements to get
  @param ix 	- indices where to get them from (in global 1d numbering) 
*/
void  PETScVector::getEntries(PetscInt ni,const PetscInt ix[], 
				      PetscScalar y[]) const
{
   VecGetValues(v, ni, ix, y);
}

/*!
    Get norm of RHS
    @param nmtype  - norm type
                     NORM_1 denotes sum_i |x_i|
                     NORM_2 denotes sqrt(sum_i (x_i)^2)
                     NORM_INFINITY denotes max_i |x_i| 
    06.2012. WW
*/
PetscReal PETScVector::getNorm(NormType  nmtype)
{
  PetscReal norm = 0.;
  VecNorm(v, nmtype, &norm); 
  return norm; 
}


void  PETScVector::restoreLocalVector(PetscScalar *loc_vec)
{
   VecRestoreArray(v, &loc_vec);
}


void PETScVector::setValue(const int i, const double value )
{

  VecSetValues(v,1,&i,&value,INSERT_VALUES);
}

void  PETScVector::setValues( PetscInt ni, const PetscInt ix[], 
                                       const PetscScalar y[],InsertMode iora) 
{
   VecSetValues(v, ni, ix, y, iora); 
}



void PETScVector::addValue(const int i, const double value,InsertMode mode )
{

  VecSetValue(v, i, value, mode);
}


void PETScVector::nullize( )
{

   VecSet(v, 0.0);
} 



void PETScVector::Viewer(std::string file_name)
{
  PetscViewer viewer;
  std::string fname = file_name + ".txt";
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, fname.c_str(), &viewer);
  PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);

  finalAssemble();


  //PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_VTK);
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
 
